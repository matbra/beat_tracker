clear, close all, clc

dir_signals = fullfile(dirup(2), 'impulse_noise', 'signals');

filename_input = 'modell_1.wav';
L_block = 1024*8;
L_feed = L_block/2;
L_DFT = L_block;
vec_window = sqrt(hann(L_block, 'periodic'));

% load the input signal
[x, fs] = wavread(fullfile(dir_signals, filename_input));
L_x = length(x);

% create the time vector
vec_t = (0:L_x-1)' / fs;

% phase transform
[mat_X, vec_f, vec_t_block] = spectrogram(x, vec_window, (L_block-L_feed), L_DFT, fs);
mat_X = ones(size(mat_X)) .* exp(j * angle(mat_X)); % set all magnitudes to one
x_phat = ispectrogram(mat_X);

L_x_phat = length(x_phat);
vec_t_phat = (0:L_x_phat-1)'/fs;

figure(1);
plot(vec_t_phat, x_phat);

% that was pre-processing...

% now perform first, coarse tempo estimation
% (with acf analysis)

T_block_acf_analysis = 4; % s
L_block_acf_analysis = floor(T_block_acf_analysis * fs);

bpm_min = 40;
bpm_max = 300;
T_delay_max = 1/bpm_min * 60;
T_delay_min = 1/bpm_max * 60;
L_delay_max = floor(T_delay_max * fs);
L_delay_min = floor(T_delay_min * fs);

vec_xAxis_bpm = 60./((L_delay_min+1 : L_delay_max)'/fs);

N_blocks = floor(L_x_phat / L_block_acf_analysis);

st_coarse_tempo_information = [];
st_refined_tempo_information = [];
vec_beats = [];

for p = 1 : N_blocks
    idx = (p-1) * L_block_acf_analysis+1 : p * L_block_acf_analysis;
    x_phat_p = x_phat(idx);
    
    %     acf = ifft(fft(abs(x_phat_p)).^2);
    acf = xcorr(abs(x_phat_p), abs(x_phat_p), L_delay_max, 'unbiased');
    
    % normalize peak at zero delay to one
    acf = acf / acf(L_delay_max + 1);
    
    % cut away irrelevant parts
    % (negative delays and delays smaller than the minimum delay)
    acf = acf(L_delay_max + 1 + L_delay_min:2*L_delay_max);
    
    % apply some smoothing
    T_smooth = 1e-2;
    L_filter_smooth = floor(T_smooth * fs);
    acf = filtfilt(1/L_filter_smooth * ones(L_filter_smooth, 1), 1, acf);
    
    % determine threshold
    th = mean(acf) + 2 * std(acf);
    
    if true
        figure(2)
        
        plot(vec_xAxis_bpm, acf);
        %         plot(acf);
        
        % draw the threshold
        line([vec_xAxis_bpm(1), vec_xAxis_bpm(end)], [th th], 'color', 'red');
        
        xlabel('bpm');
        ylabel('acf');
    end
    
    % find peaks above the threshold
    vec_b_above_threshold = acf >= th;
    
    % find connected regions
    st_regions_above_threshold = find_sections(find(vec_b_above_threshold));
    N_regions_above_threshold = length(st_regions_above_threshold);
    
    % find maxima within each region
    for b = 1 : N_regions_above_threshold
        [~, idx_max] = max(acf(st_regions_above_threshold(b).idx_start:st_regions_above_threshold(b).idx_end));
        vec_delay(b) = idx_max + L_delay_min + st_regions_above_threshold(b).idx_start-2;
    end
    
    if N_regions_above_threshold > 0
        % there seems to be some kind of tempo
        
        % store information in a struct
        st_coarse_tempo_information(end+1).idx_start = idx(1);
        st_coarse_tempo_information(end).idx_end = idx(end);
        st_coarse_tempo_information(end).vec_delay = vec_delay;
        
        % try to find the "phase" of the tempo grid
        % (acf only yields tempo, not where the measures start...)
        
        % create comb with teeth matched to the detected tempo candidates
        vec_comb_corr = zeros(length(vec_delay), 1);
        for b = 1 : length(vec_delay)
            cur_delay = vec_delay(b);
            
            T_tooth = 1e-3; % s
            L_tooth = floor(T_tooth * fs);
            L_tooth_half = floor(L_tooth/2);
            L_tooth = L_tooth_half * 2; % is now even
            
            vec_tooth = hann(L_tooth, 'symmetric');
            
            L_zeros = cur_delay - 2 * L_tooth_half;
            
            % the number of teeth that fit into the selected acf analysis
            % block length
            N_teeth = floor(L_block_acf_analysis / (L_tooth + L_zeros));
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
            th_beat_phase = mean(acf_cyclic) + 3 * std(acf_cyclic);
            vec_b_above_threshold = acf_cyclic >= th_beat_phase;
            
            if true
                % plot the correlation with the comb
                figure(3);
                
                plot(acf_cyclic);
                
                line([1 L_comb], repmat(th_beat_phase, 1, 2), 'color', 'red');
            end
            
            
            % find connected regions
            st_regions_above_threshold = find_sections(find(vec_b_above_threshold));
            
            % refine the first maximum
            [~, idx_max] = max(acf_cyclic(st_regions_above_threshold(1).idx_start:st_regions_above_threshold(1).idx_end));
            
            st_coarse_tempo_information(end).offset(b) = st_regions_above_threshold(1).idx_start + idx_max - 1;
            
            % try to refine the tempo estimation by searching for
            % transients in the vicinity of the predicted beats
            
            
            
            
            
        end
        
        % now we have collected all possible tempo candidates and their
        % respective beat offsets
        
        T_peak_search_region = 1e-2;
        L_peak_search_region = floor(T_peak_search_region * fs);
        L_peak_search_region_half = floor(L_peak_search_region / 2);
        L_peak_search_region = L_peak_search_region_half * 2 + 1; % is odd!
        
        if true
            % plot the areas where peaks are searched
            
            figure(4);
            
            plot(x_phat_p); hold on;
            
            % plot the projected peak positions
            idx_peak = st_coarse_tempo_information(p).offset(2); % watch out, only second offset!
            while idx_peak < length(x_phat_p)
                line([idx_peak idx_peak], [0 1], 'color', 'green');
                
                % draw the search area
                plot((idx_peak-L_peak_search_region_half : idx_peak + L_peak_search_region_half), ...
                    x_phat_p(idx_peak-L_peak_search_region_half : idx_peak + L_peak_search_region_half), ...
                    'color', 'red');
                
                idx_peak = idx_peak + st_coarse_tempo_information(p).vec_delay(2); % watch out, only second tempo!
            end
            
            hold off;
            
        end
        
        % increase the time resolution of the transient position by finding
        % the maximum of the hilbert envelope in the predicted vicinity of
        % the beats...
        idx_peak = st_coarse_tempo_information(p).offset(2); % watch out, only second offset!
            while idx_peak < length(x_phat_p)
                line([idx_peak idx_peak], [0 1], 'color', 'green');
                
                cur_vec_x_vicinity = x_phat_p(idx_peak-L_peak_search_region_half : idx_peak + L_peak_search_region_half);
                
                temp = hilbert(cur_vec_x_vicinity);
                cur_vec_x_vicinity_envelope = temp .* conj(temp);
                
                [~, idx_max] = max(cur_vec_x_vicinity_envelope);
                
                vec_beats(end+1) = idx(1) + idx_peak-L_peak_search_region_half + idx_max - 1;
                
                idx_peak = idx_peak + st_coarse_tempo_information(p).vec_delay(2); % watch out, only second tempo!
            end
        
        
        
    else
        % there doesn't seem to be some kind of tempo detectable
    end
    
    
end