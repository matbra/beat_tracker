function idx_beats = detect_beats_davies_standard(x, fs)

% modified from Matthew Davies' MATLAB function
% "davies_standard()" to work with vectors instead of
% loading and writing files

% download his original code
% from https://code.soundsoftware.ac.uk/projects/davies-beat-tracker/repository
% and put everything into the "davies_beat_tracker" subdirectory.

addpath('davies_beat_tracker');

% convert to mono
x = mean(x,2);

% if audio is not at 44khz resample
if fs~=44100,
  x = resample(x,44100,fs);
end

% read beat tracking parameters
p = bt_parms;

% generate the onset detection function
df = onset_detection_function(x,p);

% strip any trailing zeros
while (df(end)==0)
  df = df(1:end-1);
end

% get periodicity path
ppath = periodicity_path(df,p);

mode = 0; % use this to run normal algorithm.
% find beat locations
beats = dynamic_programming(df,p,ppath,mode);

idx_beats = round(beats * 44100);