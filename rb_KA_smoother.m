function [ smooth_est, smooth_pts ] = rb_KA_smoother( flags, params, comb_filt_pts, filt_wts_array, observs )
%RB_SMOOTHER Run Rao-Balckwellised particle smoother using Kim's
%Approximation

Ns = params.Ns;
Np = params.Np;
K = length(observs);

% Sample some filtering trajectories to initialise the smoothing particles
[~, parent] = systematic_resample(exp(filt_wts_array{end}), Ns);
smooth_pts = comb_filt_pts(parent);

% Loop through smoothing trajectories
for ss = 1:Ns
    
    fprintf(1, '*** Smoothing sequence %u.\n', ss)
    
    % Initialise sampler
    P = params.ARO;
    smooth_pts(ss).nonlin_samp(:,1:K-1) = zeros(P+1, K-1);
    
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
            
            % Linear probability
            lin_prob(jj) = 0;
            
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
        
    end
    
end

% RTS smooth to find the linear state
smooth_est = rts_particles( flags, params, smooth_pts, log(ones(Ns,1)/Ns), observs );

end

