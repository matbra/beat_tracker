clear, close all, clc

dir_signals = fullfile(dirup(2), 'impulse_noise', 'signals');
dir_signals = '/media/programme_2/remote_task_dir/common/signals/impulse_disturbance_full_until_1980';
% dir_signals = 'H:\testsignale\mirex\beattrack_train_2006\train';
% dir_signals = 'H:\testsignale\stoerungen\sammlungen_fuer_bestimme_anwendungen\signals\impulse_disturbance\withoutDisturbance';
% dir_signals = 'H:\testsignale\stoerungen\sammlungen_fuer_bestimme_anwendungen\signals\impulse_disturbance\withDisturbance';
% dir_signals = cd;
% dir_signals = 'E:\remote_task_dir\common\signals\impulse_disturbance_full';


% filename_input = 'modell_1.wav';
filename_input = 'roboter_1.wav';
% filename_input = 'patsy.wav';
filename_input = 'withoutDisturbance/02 - Association - Along Comes Mary_snippet.wav';

% filename_input = '073_-_1958_-_Don_Gibson_-_Oh,_Lonesome_Me.wav';
% filename_input = '537_-_1981_-_Ideal_-_Erschie�en.wav';
% filename_input = 'roboter_snippet.wav';
% filename_input = 'roboter_temp.wav';
% filename_input = 'metropolis_test.wav';
% filename_input = '827_-_1996_-_Bim_Sherman_-_Solid_as_a_Rock.wav';
% filename_input = 'ID8_BennyGoodmann.wav';
% filename_input = '01-AudioTrack 01.wav';
% filename_input = 'iw_clicks.wav';
% filename_input = '02-AudioTrack 02.wav';

% filename_input = 'withoutDisturbance\02 - Association - Along Comes Mary_snippet.wav';
% filename_input = 'withoutDisturbance\17 - Drafi Deutscher & His Magics - Shake Hands_snippet.wav';
% filename_input = 'withoutDisturbance\14 - Cpt. Kirk & - Geldunter_snippet.wav';
% filename_input = 'withoutDisturbance\01 - Ryan Adams - New York, New York_snippet.wav';
% filename_input = 'withoutDisturbance\02 - Chic - I Want Your Love_snippet.wav';
% filename_input = 'withDisturbance\RossReel_snippet.wav';

% filename_input = 'train8.wav';

% load the input signal
[x, fs] = wavread(fullfile(dir_signals, filename_input));
x = x(:,1);

% x = x(80*fs:end);
x = x(1:20*fs);

st_beat_detection_result = detect_beats(x, fs, [], 4);

y = bleepify(0.5*x, [st_beat_detection_result.st_beat_info.sample_pos], fs);

% plot the tempo curve
figure(1);
vec_beats = [st_beat_detection_result.st_beat_info.sample_pos];
plot(60./(diff(vec_beats)/fs));

%% plot the estimated beats
vec_t = (0:length(x)-1)'/fs;
    figure(2);
   
       
    plot(vec_t, x, 'black');
    hold on;
    for a = 1 : length(vec_beats)
        line(repmat(vec_t(vec_beats(a)), 2, 1), [-1 1], 'color', 'red', 'linewidth', 0.5);
    end