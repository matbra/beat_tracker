clear, close all, clc

dir_signals = fullfile(dirup(2), 'impulse_noise', 'signals');
% dir_signals = 'H:\testsignale\mirex\beattrack_train_2006\train';

% filename_input = 'modell_1.wav';
filename_input = 'roboter_1.wav';
% filename_input = 'patsy.wav';

% filename_input = 'train8.wav';

% load the input signal
[x, fs] = wavread(fullfile(dir_signals, filename_input));
x = x(:,1);

st_beat_detection_result = detect_beats(x, fs);

y = bleepify(x, [st_beat_detection_result.st_beat_info.sample_pos], fs);