function st_beat_detection_result = detect_beats(x, fs, T_peak_search_region, tempo_multiplier)

st_beat_detection_result.st_beat_info.sample_pos = [];

% some global parameters:
L_block = 1024*8;
L_feed = L_block/2;
L_DFT = L_block;
T_block_acf_analysis = 4; % s - size of the blocks that are used for acf analysis
bpm_min = 40; % minimum bpm
bpm_max = 600; % maximum bpm
T_smooth = 1e-2; % time constant of the (FIR) smoothing filter used for acf smoothing
if isempty(T_peak_search_region)
    T_peak_search_region = 5e-2; % size of the area around detected peaks where transients are searched in the time domain signal
end
if isempty(tempo_multiplier)
    tempo_multiplier = 1;
end
T_tooth = 1e-3; % s - width of the offset-estimation tooths (of the comb vector)
vec_window = sqrt(hann(L_block, 'periodic'));
detectorsignal_mode = 'phat';
b_prior_tempo_estimation = true;

L_tooth = floor(T_tooth * fs);
L_tooth_half = floor(L_tooth/2);
L_tooth = L_tooth_half * 2; % is now even

vec_tooth = hann(L_tooth, 'symmetric');

% use the first channel only
x = x(:,1);
L_x = length(x);

% compute detector signal
switch(detectorsignal_mode)
    case 'phat'
        Data = spectrogram(x,vec_window,(L_block-L_feed),L_DFT,fs,'yaxis');
        Data_log = 20*log10(abs(Data)+eps);
        
        % this is the phat-transform:
        x_phat = ispecgram((ones(size(Data_log)).*exp(1j*angle(Data))), L_DFT,fs);
    case 'ar'
        x_phat = ar_error(x, 16, 1024);
    case 'product_1'
        Data = spectrogram(x,vec_window,(L_block-L_feed),L_DFT,fs,'yaxis');
        Data_log = 20*log10(abs(Data)+eps);
        
        % this is the phat-transform:
        x_phat = ispecgram((ones(size(Data_log)).*exp(1j*angle(Data))), L_DFT,fs);
        
        ar_err_sig = ar_error(x, 16, 1024);
        
        x_phat = x_phat.* ar_err_sig(1:length(x_phat));
end

L_x_phat = length(x_phat);

% that was pre-processing...

% now perform first, coarse tempo estimation
% (with acf analysis)
L_block_acf_analysis = floor(T_block_acf_analysis * fs);

T_delay_max = 1/bpm_min * 60;
T_delay_min = 1/bpm_max * 60;
L_delay_max = floor(T_delay_max * fs);
L_delay_min = floor(T_delay_min * fs);

N_full_blocks = floor(L_x_phat / L_block_acf_analysis);
N_blocks = ceil(L_x_phat / L_block_acf_analysis);

% determine block sizes
% (last block is probably smaller)
vec_L_block = [L_block_acf_analysis * ones(N_full_blocks,1); L_x_phat-N_full_blocks * L_block_acf_analysis];

if vec_L_block(end) < L_delay_max
    N_blocks = N_blocks - 1;
end

st_coarse_tempo_information = [];
vec_beats = [];

idx_start = 1;

if b_prior_tempo_estimation
    % do davies beat tracking
    try
        idx_beats_davies = detect_beats_davies_standard(x, fs);
        tempo_davies = 60 ./diff(idx_beats_davies/fs);
    catch
        % error. signal was probably too short...
        idx_beats_davies = [];
        tempo_davies = [];
    end
else
    idx_beats_davies = [];
    tempo_davies = [];
end

% some parameters for the bpm pdf
sigma_bpm = 10;

for p = 1 : N_blocks
    idx = idx_start : idx_start + vec_L_block(p) - 1;
    idx_start = idx_start + vec_L_block(p);
    x_phat_p = x_phat(idx);
    
    % skip all calculation of rhythmness for the last block
    % -> assume that the last tempo still holds...
    
    if p < N_blocks
        acf = xcorr(abs(x_phat_p), abs(x_phat_p), L_delay_max, 'unbiased');
        
        % normalize peak at zero delay to one
        acf = acf / acf(L_delay_max + 1);
        
        % cut away irrelevant parts
        % (negative delays and delays smaller than the minimum delay)
        acf = acf(L_delay_max + 1 + L_delay_min:2*L_delay_max);
        
        % apply some smoothing
        L_filter_smooth = floor(T_smooth * fs);
        acf = filtfilt(1/L_filter_smooth * ones(L_filter_smooth, 1), 1, acf);
        
        % determine threshold
        th = mean(acf) + 1 * std(acf);
                
        % find peaks above the threshold
        vec_b_above_threshold = acf >= th;
        
        % find connected regions
        st_regions_above_threshold = find_sections(find(vec_b_above_threshold));
        N_regions_above_threshold = length(st_regions_above_threshold);
        
        % find maxima within each region
        vec_delay = [];
        for b = 1 : N_regions_above_threshold
            [~, idx_max] = max(acf(st_regions_above_threshold(b).idx_start:st_regions_above_threshold(b).idx_end));
            vec_delay(b) = idx_max + L_delay_min + st_regions_above_threshold(b).idx_start-2;
        end
        
    end
    
    % apply tempo multiplicator
    % (to include sporadic rolls etc.)
    vec_delay = round(vec_delay / tempo_multiplier);
    
    % now use the davies beat tracker information to find most probable
    % tempo (delay) candidate
    if b_prior_tempo_estimation
        % determine which beats (of davies output) are within the current block
        vec_b_beats_davies_within_cur_block = idx_beats_davies >= idx(1) & ...
            idx_beats_davies <= idx(end);
        
        idx_beats_davies_p = idx_beats_davies(vec_b_beats_davies_within_cur_block);
        
        tempo_davies_p = 60 ./diff(idx_beats_davies_p/fs);
        
        tempo_davies_p_mean = median(tempo_davies_p);
        
        % generate the apriori tempo pdf
        cur_bpm_candidate = tempo_davies_p_mean/1;
        vec_pdf_bpm = zeros(bpm_max-bpm_min+1, 1);
        vec_bpm = (bpm_min:bpm_max)';
        while cur_bpm_candidate < bpm_max
            vec_pdf_bpm = vec_pdf_bpm + ...
                1/sqrt(2*pi*sigma_bpm^2) * exp(-(vec_bpm - cur_bpm_candidate).^2/(2*sigma_bpm^2));
            cur_bpm_candidate = cur_bpm_candidate * 2; % assuming only doubling-errors (no triplets etc.)
        end
        
        vec_pdf_bpm = normalize_pdf(vec_pdf_bpm, vec_bpm);
        
        % determine the probability of all tempo candidates
        vec_p_bpm = [];
        for a = 1 : N_regions_above_threshold
            cur_tempo_candidate_from_delay = 60/(vec_delay(a) /fs);
            vec_p_bpm(a) = interp1(vec_bpm, vec_pdf_bpm, cur_tempo_candidate_from_delay, 'linear');
        end
        
        % threshold the tempo candidates
        th_p_bpm = 0.006;
        vec_b_probable_tempo = vec_p_bpm > th_p_bpm;
        
        % only allow the fastest tempo
        temp_idx_fastest = find(vec_b_probable_tempo, 1, 'first'); % TODO: shouldn't this be 'last' ??
        vec_b_probable_tempo = false(N_regions_above_threshold, 1);
        vec_b_probable_tempo(temp_idx_fastest) = true;
    else
        temp_idx_fastest = 1;%N_regions_above_threshold; % TODO: is it correct to pick the last one?
        vec_b_probable_tempo = true(N_regions_above_threshold, 1);
    end
    
    st_coarse_tempo_information(end+1).b_valid = false;
    
    if N_regions_above_threshold > 0 && any(vec_b_probable_tempo)
        % there seems to be some kind of tempo
        
        % store information in a struct
        st_coarse_tempo_information(end).idx_start = idx(1);
        st_coarse_tempo_information(end).idx_end = idx(end);
        st_coarse_tempo_information(end).vec_delay = vec_delay;
        st_coarse_tempo_information(end).b_valid = true;
        
        % try to find the "phase" of the tempo grid
        % (acf only yields tempo, not where the measures start...)
        
        % create comb with teeth matched to the detected tempo candidates
        for b = 1 : length(vec_delay)
            
            if ~vec_b_probable_tempo(b), continue; end
            
            cur_delay = vec_delay(b);
            
            L_zeros = cur_delay - 2 * L_tooth_half;
            
            % the number of teeth that fit into the selected acf analysis
            % block length
            N_teeth = floor(vec_L_block(p) / (L_tooth + L_zeros));
            L_comb = N_teeth * (L_tooth + L_zeros);
            
            % generate the comb
            vec_comb = zeros(L_comb, 1);
            for c = 1 : N_teeth
                vec_comb((c-1) * (L_tooth + L_zeros) + 1 : (c-1) * (L_tooth + L_zeros) + L_tooth) = vec_tooth;
            end
            
            acf_cyclic = ifft(fft(abs(x_phat_p(1:L_comb))) .* fft(vec_comb));
            
            % smooth...
            acf_cyclic = filtfilt(1/L_filter_smooth * ones(L_filter_smooth ,1), 1, acf_cyclic);
            
            % find the first maximum to determine the offset for this tempo
            % candidate
            fac_std = 2;
            vec_b_above_threshold = false(L_comb, 1);
            while ~nnz(vec_b_above_threshold)
                th_beat_phase = mean(acf_cyclic) + fac_std * std(acf_cyclic);
                vec_b_above_threshold = acf_cyclic >= th_beat_phase;
                fac_std = fac_std * 0.9;
                
                if isnan(th_beat_phase)
                    break;
                end
            end
            
            % find connected regions
            st_regions_above_threshold = find_sections(find(vec_b_above_threshold));
            
            % refine the first maximum
            if ~isempty(st_regions_above_threshold)
                [~, idx_max] = max(acf_cyclic(st_regions_above_threshold(1).idx_start:st_regions_above_threshold(1).idx_end));
                st_coarse_tempo_information(end).offset(b) = st_regions_above_threshold(1).idx_start + idx_max - 1;
            else
                idx_max = 1; % todo: this is evil!
                st_coarse_tempo_information(end).offset(b) = 0;
            end
            
            % try to refine the tempo estimation by searching for
            % transients in the vicinity of the predicted beats
        end
        
        % now we have collected all possible tempo candidates and their
        % respective beat offsets
        
        L_peak_search_region = floor(T_peak_search_region * fs);
        L_peak_search_region_half = floor(L_peak_search_region / 2);
        L_peak_search_region = L_peak_search_region_half * 2 + 1; % is odd!
        
        % increase the time resolution of the transient position by finding
        % the maximum of the hilbert envelope in the predicted vicinity of
        % the beats...
        idx_peak = st_coarse_tempo_information(p).offset(temp_idx_fastest); % watch out, only second offset!
        while idx_peak < length(x_phat_p)
            
            idx_start_search = idx_peak-L_peak_search_region_half;
            
            if idx_start_search < 1
                idx_start_search = 1;
            end
            
            idx_end_search = idx_peak + L_peak_search_region_half;
            if idx_end_search > vec_L_block(p)
                idx_end_search = vec_L_block(p);
            end
            
            cur_vec_x_vicinity = x_phat_p(idx_start_search : idx_end_search);
            
            temp = hilbert(cur_vec_x_vicinity);
            cur_vec_x_vicinity_envelope = temp .* conj(temp);
            
            [~, idx_max] = max(cur_vec_x_vicinity_envelope);
            
            vec_beats(end+1) = idx(1) + idx_start_search + idx_max - 1;
            
            idx_peak = idx_peak + st_coarse_tempo_information(p).vec_delay(temp_idx_fastest); % watch out, only second tempo!
        end
    else
        % there doesn't seem to be some kind of tempo detectable
    end
end

% prepare the return struct
for a = 1 : length(vec_beats)
    st_beat_detection_result.st_beat_info(a).beat_index = a;
    st_beat_detection_result.st_beat_info(a).sample_pos = vec_beats(a);
end

% some global information
st_beat_detection_result.st_global_info.idx_start_analyse = 1;
st_beat_detection_result.st_global_info.idx_end_analyse = sum(vec_L_block);%N_blocks * L_block_acf_analysis;