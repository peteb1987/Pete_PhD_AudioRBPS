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
flags.data = 2;
params.procvar_decay = 0.995;           % Decay coefficient on log of process variance
params.logprocvar_vr = 0.1;5E-3;            % Transition variance of log of process variance
params.ref_trans_vr = 5E-3;             % Transition variance of reflection coefficients
params.noise_vr = 0.001;                % Noise variance
params.min_log_proc_var = -10;          % Lower limit on log(noise variance)
params.ARO = 5;                         % Order of AR model
params.Np = 100;                        % Number of filtering particles
params.Ns = 10;                         % Number of smoothing trajectories
params.resam_thresh = 1;              % Proportion of Np which ESS must exceed for resampling to occur
params.init_ref_vr = 0.5;               % Prior variance for reflection coefficients
params.init_logprocvar_mn = -8;         % Prior mean for log of process variance
params.init_logprocvar_vr = 0.1;        % Prior variance for log of process variance

params.K = 500;                        % Number of samples
params.fs = 44100;                      % Input audio sampling frequency
params.ds = 4;                          % Downsampling rate

if flags.data == 1
    
    % Get an audio wav file
    [true_audio] = wavread(filename);
    
    % Downsample
    ds = params.ds;
    B = fir1(100, 1/ds);
    true_audio = filter(B, 1, true_audio);
    true_audio = true_audio(1:ds:end);
    params.fs = params.fs/ds;
    
    % Truncate
    start_samp = round(params.fs/2)+250;
    num_samp = params.K;
    true_audio = true_audio(start_samp+1:start_samp+num_samp);
    
elseif flags.data == 2
    
    params.fs = params.fs/params.ds;
    
    % Generate some data
    [true_audio, true_ar_coeffs, true_proc_var] = generate_data(flags, params);
    
end

% Add some noise
noisy_audio = true_audio + mvnrnd(zeros(size(true_audio)), params.noise_vr);

%% Run filter
fprintf(1, 'Filtering\n')
[ init_pts ] = initialise_particles(flags, params);
[ filt_est, comb_filt_pts, final_filt_pts, filt_wts_array ] = rb_filter( flags, params, init_pts, noisy_audio );

%% RTS for Kitagawa smoothed estimate
fprintf(1, 'RTS smoothing\n')
[ kita_est ] = rts_particles( flags, params, final_filt_pts, filt_wts_array{end}, noisy_audio );

%% Run smoother
fprintf(1, 'New-type smoothing\n')
[ smooth_est, smooth_pts ] = rb_smoother( flags, params, comb_filt_pts, filt_wts_array, noisy_audio );

%% Run Fong-type sampling smoother
fprintf(1, 'Fong-type smoothing\n')
[ samp_smooth_est, samp_smooth_pts ] = rb_sampling_smoother( flags, params, comb_filt_pts, filt_wts_array, noisy_audio );

%% Run Kim's approximation sampling smoother
fprintf(1, 'Kim''s approximation smoothing\n')
[ KA_smooth_est, KA_smooth_pts ] = rb_KA_smoother( flags, params, comb_filt_pts, filt_wts_array, noisy_audio );

%% Measure SNRs
input_SNR = SNR(true_audio, noisy_audio)
filt_SNR = SNR(true_audio, filt_est)
kita_SNR = SNR(true_audio, kita_est)
smooth_SNR = SNR(true_audio, smooth_est)
samp_smooth_SNR = SNR(true_audio, samp_smooth_est)
KA_smooth_SNR = SNR(true_audio, KA_smooth_est)

%% Play
player = audioplayer(true_audio, params.fs);
play(player);% pause;
player = audioplayer(noisy_audio, params.fs);
play(player);% pause;
player = audioplayer(filt_est, params.fs);
play(player);% pause;
player = audioplayer(kita_est, params.fs);
play(player);% pause;
player = audioplayer(smooth_est, params.fs);
play(player);% pause;
player = audioplayer(samp_smooth_est, params.fs);
play(player);% pause;
player = audioplayer(KA_smooth_est, params.fs);
play(player);% pause;

%% Plotting
figure, hold on, plot(true_audio, 'b'); plot(noisy_audio, 'r');
figure, hold on, plot(true_audio, 'b'); plot(filt_est, 'm');
figure, hold on, plot(true_audio, 'b'); plot(kita_est, 'c');
figure, hold on, plot(true_audio, 'b'); plot(smooth_est, 'g');
figure, hold on, plot(true_audio, 'b'); plot(samp_smooth_est, 'k');
figure, hold on, plot(true_audio, 'b'); plot(KA_smooth_est, 'color', [.5 .5 .5]);

%% Nonlinear parts
filt_pts_nonlin_samps = cat(3, comb_filt_pts.nonlin_samp);
figure, hold on, plot(log(squeeze(filt_pts_nonlin_samps(end,:,:))));
if exist('true_proc_var','var')
    plot(log(true_proc_var), 'r', 'linewidth', 2);
end
figure, hold on, plot(squeeze(filt_pts_nonlin_samps(1,:,:)));
if exist('true_ar_coeffs','var')
    plot(true_ar_coeffs(1,:), 'r', 'linewidth', 2);
end

final_filt_pts_nonlin_samps = cat(3, final_filt_pts.nonlin_samp);
figure, hold on, plot(log(squeeze(final_filt_pts_nonlin_samps(end,:,:))));
if exist('true_proc_var','var')
    plot(log(true_proc_var), 'r', 'linewidth', 2);
end
figure, hold on, plot(squeeze(final_filt_pts_nonlin_samps(1,:,:)));
if exist('true_ar_coeffs','var')
    plot(true_ar_coeffs(1,:), 'r', 'linewidth', 2);
end

smooth_pts_nonlin_samps = cat(3, smooth_pts.nonlin_samp);
figure, hold on, plot(log(squeeze(smooth_pts_nonlin_samps(end,:,:))));
if exist('true_proc_var','var')
    plot(log(true_proc_var), 'r', 'linewidth', 2);
end
figure, hold on, plot(squeeze(smooth_pts_nonlin_samps(1,:,:)));
if exist('true_ar_coeffs','var')
    plot(true_ar_coeffs(1,:), 'r', 'linewidth', 2);
end

samp_smooth_pts_nonlin_samps = cat(3, samp_smooth_pts.nonlin_samp);
figure, hold on, plot(log(squeeze(samp_smooth_pts_nonlin_samps(end,:,:))));
if exist('true_proc_var','var')
    plot(log(true_proc_var), 'r', 'linewidth', 2);
end
figure, hold on, plot(squeeze(samp_smooth_pts_nonlin_samps(1,:,:)));
if exist('true_ar_coeffs','var')
    plot(true_ar_coeffs(1,:), 'r', 'linewidth', 2);
end

KA_smooth_pts_nonlin_samps = cat(3, KA_smooth_pts.nonlin_samp);
figure, hold on, plot(log(squeeze(KA_smooth_pts_nonlin_samps(end,:,:))));
if exist('true_proc_var','var')
    plot(log(true_proc_var), 'r', 'linewidth', 2);
end
figure, hold on, plot(squeeze(KA_smooth_pts_nonlin_samps(1,:,:)));
if exist('true_ar_coeffs','var')
    plot(true_ar_coeffs(1,:), 'r', 'linewidth', 2);
end

%% Measure nonlinear error

if exist('true_ar_coeffs','var')
    filt_mean_nonlin = mean(cat(3, final_filt_pts.nonlin_samp), 3);
    filt_ar_rmse = sqrt( mean(mean((filt_mean_nonlin(1:params.ARO,:)-true_ar_coeffs).^2)) )
    filt_pv_rmse = sqrt( mean((filt_mean_nonlin(params.ARO+1,:)-true_proc_var).^2) )
    
    smooth_mean_nonlin = mean(cat(3, smooth_pts.nonlin_samp), 3);
    smooth_ar_rmse = sqrt( mean(mean((smooth_mean_nonlin(1:params.ARO,:)-true_ar_coeffs).^2)) )
    smooth_pv_rmse = sqrt( mean((smooth_mean_nonlin(params.ARO+1,:)-true_proc_var).^2) )
    
    samp_smooth_mean_nonlin = mean(cat(3, samp_smooth_pts.nonlin_samp), 3);
    samp_smooth_ar_rmse = sqrt( mean(mean((samp_smooth_mean_nonlin(1:params.ARO,:)-true_ar_coeffs).^2)) )
    samp_smooth_pv_rmse = sqrt( mean((samp_smooth_mean_nonlin(params.ARO+1,:)-true_proc_var).^2) )
    
    KA_smooth_mean_nonlin = mean(cat(3, KA_smooth_pts.nonlin_samp), 3);
    KA_smooth_ar_rmse = sqrt( mean(mean((KA_smooth_mean_nonlin(1:params.ARO,:)-true_ar_coeffs).^2)) )
    KA_smooth_pv_rmse = sqrt( mean((KA_smooth_mean_nonlin(params.ARO+1,:)-true_proc_var).^2) )
end

%% Count unique particles
[filt_Nup, filt_Nut] = count_unique_particles(params.K, final_filt_pts);
[smooth_Nup, smooth_Nut] = count_unique_particles(params.K, smooth_pts);
[samp_smooth_Nup, samp_smooth_Nut] = count_unique_particles(params.K, samp_smooth_pts);
[KA_smooth_Nup, KA_smooth_Nut] = count_unique_particles(params.K, KA_smooth_pts);
mean(filt_Nup)
mean(smooth_Nup)
mean(samp_smooth_Nup)
mean(KA_smooth_Nup)
figure, hold on, plot(filt_Nup, 'm'), plot(smooth_Nup, 'g'), plot(samp_smooth_Nup, 'k'), plot(KA_smooth_Nup, 'color', [.5 .5 .5])