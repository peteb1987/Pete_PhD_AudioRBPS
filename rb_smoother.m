function [ smooth_est, smooth_pts ] = rb_smoother( flags, params, comb_filt_pts, filt_wts_array, observs )
%RB_SMOOTHER Run Rao-Balckwellised particle smoother

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
    
    % Delete everything before the final frame
    P = params.ARO;
    smooth_pts(ss).nonlin_samp(:,1:K-P) = zeros(P+1, K-P);
    smooth_pts(ss).lin_mn(:,1:K-P) = zeros(P, K-P);
    smooth_pts(ss).lin_vr(:,:,1:K-P) = zeros(P, P, K-P);
    
    % MAKE THIS INITIALISATION BETTER!!!!!
    
    % Initialise backward filter
    back_mn = smooth_pts(ss).lin_mn(:,K);
    back_vr = smooth_pts(ss).lin_vr(:,:,K);
    
    % Loop backwards in time
    for kk = K-params.ARO:-1:1
        
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

            % P-step forward prediction
            forw_mn = comb_filt_pts(jj).lin_mn(1:P,kk);
            forw_vr = comb_filt_pts(jj).lin_vr(1:P,1:P,kk);
            for pp = 1:P
                [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+pp), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+pp));
                [forw_mn, forw_vr] = kf_predict(forw_mn, forw_vr, A, Q);
            end
            forw_vr = (forw_vr+forw_vr')/2;
            
            % Linear probability
            test_vr = back_vr+forw_vr;
            test_mn = back_mn - forw_mn;
            lin_prob(jj) = log( mvnpdf( zeros(1,P), test_mn', test_vr) );
            
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
        
        % Backwards Kalman smoother update
        [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+P), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+P));
        [back_mn, back_vr] = kf_predict(back_mn, back_vr, inv(A), A\Q/A');
        H = [zeros(1, P-1) 1]; R = params.noise_vr;
        [back_mn, back_vr] = kf_update_1D(back_mn, back_vr, observs(kk), H, R);
        back_vr = (back_vr+back_vr')/2;
        
        assert(isposdef(back_vr));
        
    end
    
end

% RTS smooth to find the linear state
smooth_est = rts_particles( flags, params, smooth_pts, log(ones(Ns,1)/Ns), observs );

end

