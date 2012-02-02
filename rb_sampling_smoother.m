function [ smooth_est, smooth_pts ] = rb_sampling_smoother( flags, params, comb_filt_pts, filt_wts_array, observs )
%RB_SMOOTHER Run Rao-Balckwellised particle smoother of Fong 2004, with
%linear sampling

Ns = params.Ns;
Np = params.Np;
K = length(observs);
est = zeros(size(observs));

% Sample some filtering trajectories to initialise the smoothing particles
[~, parent] = systematic_resample(exp(filt_wts_array{end}), Ns);
smooth_pts = comb_filt_pts(parent);

% Loop through smoothing trajectories
for ss = 1:Ns
    
    fprintf(1, '*** Smoothing sequence %u.\n', ss)
    
    % Initialise sampler
    P = params.ARO;
    smooth_pts(ss).nonlin_samp(:,1:K-1) = zeros(P+1, K-1);
    smooth_pts(ss).lin_samp = zeros(P, K-1);
    smooth_pts(ss).lin_samp(:,K) = mvnrnd( smooth_pts(ss).lin_mn(:,K)', smooth_pts(ss).lin_vr(:,:,K) )';
    
    % Initialise backward filter
    back_mn = smooth_pts(ss).lin_mn(:,K);
    back_vr = smooth_pts(ss).lin_vr(:,:,K);
    
    % Loop backwards in time
    for kk = K-1:-1:1
        
        if mod(kk,100)==0
            fprintf(1, '*** Time point %u.\n', kk)
        end
        
%         P = min(params.ARO, kk);
        
        % Calculate the sampling weights
        sampling_weights = zeros(Np, 1);
        trans_prob = zeros(Np, 1);
        lin_prob = zeros(Np, 1);
        
        % Loop through the previous frame particles
        for jj = 1:Np
            
            % Calculate nonlinear transition density
            [~, trans_prob(jj)] = sample_nonlin_transdens(flags, params, comb_filt_pts(jj).nonlin_samp(:,kk), smooth_pts(ss).nonlin_samp(:,kk+1));

            % Forward prediction
            forw_mn = comb_filt_pts(jj).lin_mn(1:P,kk);
            forw_vr = comb_filt_pts(jj).lin_vr(1:P,1:P,kk);
            [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+1), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+1));
            [pred_forw_mn, pred_forw_vr] = kf_predict(forw_mn, forw_vr, A, Q);
            pred_forw_vr = (pred_forw_vr+pred_forw_vr')/2;
            
            % Linear probability
            lin_prob(jj) = log( mvnpdf( smooth_pts(ss).lin_samp(:,kk+1)', pred_forw_mn', pred_forw_vr) );
            
            % Add up weights
            sampling_weights(jj) = filt_wts_array{kk}(jj) ...
                + trans_prob(jj) ...
                + lin_prob(jj);
            
        end
        
        % Normalise weights
        sampling_weights = sampling_weights - max(sampling_weights);
        lin_weights = exp(sampling_weights); lin_weights = lin_weights/sum(lin_weights);
        sampling_weights = log(lin_weights);
        
        % Backwards sampling
        samp_ind = randsample(Np, 1, true, exp(sampling_weights));
        smooth_pts(ss).nonlin_samp(:, kk) = comb_filt_pts(samp_ind).nonlin_samp(:, kk);
        
        % Kalman smoother
        [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+1), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+1));
        smooth_gain = forw_vr * A * pred_forw_vr;
        back_mn = forw_mn + smooth_gain*(back_mn - pred_forw_mn);
        back_vr = forw_vr + smooth_gain*(back_vr - pred_forw_vr)*smooth_gain';
        back_vr = (back_vr+back_vr')/2;
        
        assert(isposdef(back_vr));
        
        % Linear sampling
        smooth_pts(ss).lin_samp(:, kk) = mvnrnd(back_mn', back_vr)';
        
    end
    
end

% RTS smooth to find the linear state
smooth_est = rts_particles( flags, params, smooth_pts, log(ones(Ns,1)/Ns), observs );

end

