function [ audio_output, filt_pts, smooth_pts ] = audio_noise_reduction( flags, params, audio_input, true_audio )
%AUDIO_SMOOTHER - Run noise removal smoothing algorithm

Nsamp = length(audio_input);

% Create some initial particles
init_ref = inf(params.ARO,params.Np);
for ii = 1:params.Np
for pp = 1:params.ARO
    while abs(init_ref(pp,ii))>1
        init_ref(pp,ii) = normrnd(0, sqrt(params.init_ref_vr));
    end
end
end
init_ar = num2cell([ step_up(init_ref);
                     exp(normrnd(params.init_logprocvar_mn*ones(1,params.Np), sqrt(params.init_logprocvar_vr)))] , 1)';
filt_pts = struct('nonlin_samp', init_ar, 'lin_mn', zeros(params.ARO,1), 'lin_vr', 1E-20*eye(params.ARO));

% Block loop
NB = max(ceil(Nsamp/params.block_length), 1);
filt_pts_arr = cell(NB,1);
for bb = 1%:NB
    
    bb
    
    % Work out where block starts and finishes
    block_start = (bb-1)*params.block_length+1;
    block_end = min( bb*params.block_length, Nsamp );
    
    last_state_prev_pts = filt_pts;
    for ii = 1:params.Np
        last_state_prev_pts(ii).nonlin_samp(:,1:end-1) = [];
        last_state_prev_pts(ii).lin_mn(:,1:end-1) = [];
        last_state_prev_pts(ii).lin_vr(:,:,1:end-1) = [];
    end
    
    % Run filter
    [ filt_est, final_filt_pts, filt_pts, filt_weights ] = rb_filter( flags, params, last_state_prev_pts, audio_input(block_start:block_end) );
    [ final_filt_pts ] = rts_particles(flags, params, final_filt_pts);
    
    filt_SNR = SNR(true_audio, filt_est);
    figure, hold on, plot(true_audio), plot(filt_est, 'r');
    filt_player = audioplayer(filt_est, params.fs);
    play(filt_player);
    
    KitSmooth_SNR = SNR(true_audio, final_filt_pts(1).smooth_mn(1,:)');
    figure, hold on, plot(true_audio), plot(final_filt_pts(1).smooth_mn(1,:), 'r');
    KitSmooth_player = audioplayer(final_filt_pts(1).smooth_mn(1,:), params.fs);
    play(KitSmooth_player);
    
    nonlin_samps = mean(cat(3, final_filt_pts.nonlin_samp), 3);
    figure, plot(nonlin_samps(1:end-1,:)')
    figure, plot(log(nonlin_samps(end,:)))
    
    filt_pts_arr{bb} = final_filt_pts;
    
    %Run smoother
    [ smooth_pts ] = rb_smoother( flags, params, filt_pts, filt_weights, audio_input(block_start:block_end) );
    [ smooth_pts ] = rts_particles(flags, params, smooth_pts, audio_input, zeros(params.ARO,1), 1E-20*eye(params.ARO));
    
    smooth_lin_mn = mean(cat(3, smooth_pts.smooth_mn), 3);
    RBPS_SNR = SNR(true_audio, smooth_lin_mn(1,:)');
    figure, hold on, plot(true_audio), plot(smooth_lin_mn(1,:), 'r');
    RBPS_player = audioplayer(smooth_lin_mn(1,:), params.fs);
    play(RBPS_player);
    
    % Join block smoothing trajectories
    
    % Throw away intermediate particles
    
    
end

% lin_pts = cat(3, filt_pts.lin_mn);
% best_est = squeeze(mean(lin_pts,3));
% nonlin_pts = cat(3, filt_pts.nonlin_samp);
% audio_output = [best_est(end, params.ARO:end), flipud(best_est(1:params.ARO-1, end))']';
% audio_output = mean(cat(3, smooth_pts.back_lin_mn), 3);
% audio_output(2:end,:) = [];

audio_output = 0;

end

