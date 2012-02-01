% Rao-Blackwellised particle smoothing for audio de-noising

clup
dbstop if error

%% Set-up

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

% Set inputs
filename = 'TH_TED.wav';

% Set parameters
flags = [];
params.procvar_decay = 0.995;           % Decay coefficient on log of process variance
params.logprocvar_vr = 5E-3;            % Transition variance of log of process variance
params.ref_trans_vr = 5E-3;             % Transition variance of reflection coefficients
params.noise_vr = 0.0002;               % Noise variance
params.ARO = 6;                         % Order of AR model
params.Np = 100;                        % Number of filtering particles
params.Ns = 10;                         % Number of smoothing trajectories
params.resam_thresh = 0.5;              % Proportion of Np which ESS must exceed for resampling to occur
params.init_ref_vr = 0.5;               % Prior variance for reflection coefficients
params.init_logprocvar_mn = -5;         % Prior mean for log of process variance
params.init_logprocvar_vr = 0.1;        % Prior variance for log of process variance

params.K = 1000;                        % Number of samples
params.fs = 44100;                      % Input audio sampling frequency
params.ds = 4;                          % Downsampling rate

% Get an audio wav file
[true_audio] = wavread(filename);

% Downsample
ds = params.ds;
B = fir1(100, 1/ds);
true_audio = filter(B, 1, true_audio);
true_audio = true_audio(1:ds:end);
params.fs = params.fs/ds;

% Truncate
start_samp = round(params.fs/2);
num_samp = params.K;
true_audio = true_audio(start_samp+1:start_samp+num_samp);

% Add some noise
noisy_audio = true_audio + mvnrnd(zeros(size(true_audio)), params.noise_vr);

%% Run filter
[ init_pts ] = initialise_particles(flags, params);
[ filt_est, comb_filt_pts, final_filt_pts, filt_wts_array ] = rb_filter( flags, params, init_pts, noisy_audio );

%% RTS for Kitagawa smoothed estimate
[ kita_est ] = rts_particles( flags, params, final_filt_pts, filt_wts_array{end}, noisy_audio );

%% Run smoother
[ smooth_est, smooth_pts ] = rb_smoother( flags, params, comb_filt_pts, filt_wts_array, noisy_audio );

%% Measure SNRs
input_SNR = SNR(true_audio, noisy_audio);
filt_SNR = SNR(true_audio, filt_est);
kita_SNR = SNR(true_audio, kita_est);
smooth_SNR = SNR(true_audio, smooth_est);

%% Play
player = audioplayer(true_audio, params.fs);
play(player)
player = audioplayer(noisy_audio, params.fs);
play(player)
player = audioplayer(filt_est, params.fs);
play(player)
player = audioplayer(kita_est, params.fs);
play(player)
player = audioplayer(smooth_est, params.fs);
play(player)

%% Plotting
figure, hold on, plot(true_audio, 'b'); plot(noisy_audio, 'r');
figure, hold on, plot(true_audio, 'b'); plot(filt_est, 'm');
figure, hold on, plot(true_audio, 'b'); plot(kita_est, 'c');
figure, hold on, plot(true_audio, 'b'); plot(smooth_est, 'g');

%% Nonlinear parts
filt_pts_nonlin_samps = cat(3, comb_filt_pts.nonlin_samp);
figure, plot(log(squeeze(filt_pts_nonlin_samps(end,:,:))));
figure, plot(squeeze(filt_pts_nonlin_samps(1,:,:)));

smooth_pts_nonlin_samps = cat(3, smooth_pts.nonlin_samp);
figure, plot(log(squeeze(smooth_pts_nonlin_samps(end,:,:))));
figure, plot(squeeze(smooth_pts_nonlin_samps(1,:,:)));
