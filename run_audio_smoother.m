clup
dbstop if error

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

% Set inputs
filename = 'TH_TED.wav';

% Set parameters
flags.test = true;
params.fs = 44100;                      % Audio sampling frequency
params.noise_vr = 0.0002;               % Noise variance
params.ARO = 1;                         % Order of AR model
params.Np = 20;                         % Number of filtering particles
params.Ns = 10;                         % Number of smoothing trajectories
params.block_length = 1E10;              % Length of processing blocks
params.init_ref_vr = 0.5;              % Prior variance for reflection coefficients
params.init_logprocvar_mn = -5;       % Prior mean for log of process variance
params.init_logprocvar_vr = 0.1;       % Prior variance for log of process variance
params.procvar_decay = 0.995;          % Decay coefficient on log of process variance
params.logprocvar_vr = 5E-3;            % Transition variance of log of process variance
params.ref_trans_vr = 5E-3;            % Transition variance of reflection coefficients
params.K = ceil(0.05*params.fs);

% Get an audio wav file
[true_audio] = wavread(filename);

% Downsample
ds = 4;
B = fir1(100, 0.25);
true_audio = filter(B, 1, true_audio);
true_audio = true_audio(1:ds:end);
params.fs = params.fs/ds;
params.K = ceil(params.K/ds);

% Truncate
offset = round(params.fs/2);
true_audio = true_audio(offset+1:offset+params.K);

%%

% Add some noise
noisy_audio = true_audio + mvnrnd(zeros(size(true_audio)), params.noise_vr);

% figure, plot(true_audio(6*params.block_length+1:7*params.block_length))
% figure, plot(noisy_audio(6*params.block_length+1:7*params.block_length))

player = audioplayer(true_audio, params.fs);
play(player)

% Measure SNR
input_SNR = SNR(true_audio, noisy_audio);

%%

% Call filter/smoother algorithm
[est_audio, filt_pts, smooth_pts] = audio_noise_reduction(flags, params, noisy_audio, true_audio );

% % Measure SNR
% output_SNR = SNR(true_audio, est_audio);

%%
output_SNR = SNR(true_audio, est_audio);
input_SNR = SNR(true_audio, noisy_audio);

%%
% figure, hold on, plot(true_audio, 'b'), plot(noisy_audio, 'r'), plot(est_audio, 'g')
% nonlin_pts = mean(cat(3, pts.nonlin_samp), 3);
% figure, plot(log(nonlin_pts(end,:))');

%%
% est_player = audioplayer(est_audio, params.fs);
% noisy_player = audioplayer(noisy_audio, params.fs);
%%
% play(noisy_player)
%%
% play(est_player)