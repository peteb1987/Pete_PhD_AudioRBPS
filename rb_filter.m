function [ est, pts_store, final_filt_pts, wts_array ] = rb_filter( flags, params, init_pts, observs )
%RB_FILTER Rao-Blackwellised particle filter

% Initialise some things
Np = params.Np;
K = length(observs);

pts = init_pts;
wts_array = cell(K,1);
est = zeros(K,1);

preESS = zeros(K,1);
postESS = zeros(K,1);
num_resam = 0;
last_weights = zeros(Np, 1);

% Extend particles
for ii = 1:Np
    pts(ii).lin_vr = cat(3, pts(ii).lin_vr, repmat(1E-20*eye(params.ARO), [1 1 K-1]));
    pts(ii).lin_mn = [pts(ii).lin_mn, repmat(zeros(params.ARO,1), [1 K-1])];
    pts(ii).nonlin_samp = [pts(ii).nonlin_samp, repmat(zeros(params.ARO+1,1), [1 K-1])];
end
pts_store = pts;

tic;

% Loop through time
for kk = 1:K
    
    if mod(kk,100)==0
        fprintf(1, '*** Time point %u. Last 100 took %f s with %u resampling steps.\n', kk, toc, num_resam)
        num_resam = 0;
        tic;
    end
    
    weights = zeros(Np, 1);
    
    % Get the index of the last time. This ensures that the time 0 prior is overwritten
    last_idx = max(1,kk-1);
    P = min(params.ARO, kk);
    
    % Loop through particles
    for ii = 1:Np
        
        prev_nonlin_samp = pts(ii).nonlin_samp(:,last_idx);
        prev_lin_mn = pts(ii).lin_mn(1:P,last_idx);
        prev_lin_vr = pts(ii).lin_vr(1:P,1:P,last_idx);
        
        % Propose new nonlinear state (from transition density)
        pts(ii).nonlin_samp(:,kk) = sample_nonlin_transdens(flags, params, prev_nonlin_samp);
        
        % Run Kalman filter
        [A, Q] = construct_transmats(flags, params, pts(ii).nonlin_samp(1:P,kk), pts(ii).nonlin_samp(params.ARO+1,kk));
        [pred_mn, pred_vr] = kf_predict(prev_lin_mn, prev_lin_vr, A, Q);
        
        H = [1 zeros(1, P-1)];
        R = params.noise_vr;
        [pts(ii).lin_mn(1:P,kk), pts(ii).lin_vr(1:P,1:P,kk), ~, ~, ~, pred_lhood] = kf_update_1D(pred_mn, pred_vr, observs(kk), H, R);
        pts(ii).lin_vr(1:P,1:P,kk) = (pts(ii).lin_vr(1:P,1:P,kk)+pts(ii).lin_vr(1:P,1:P,kk)')/2;
        
        assert(isposdef(pred_vr));
        assert(isreal(pred_lhood));
        assert(isposdef(pts(ii).lin_vr(1:P,1:P,kk)));
        
        % Update weight
        weights(ii) = last_weights(ii) + log(pred_lhood);
        
    end
    
    % Normalise weights
    lin_weights = exp(weights); lin_weights = lin_weights/sum(lin_weights);
    weights = log(lin_weights);
    
    % Systematic resampling
    preESS(kk) = ESS(weights);
    if preESS(kk)<Np * params.resam_thresh
        [ ~, parent ] = systematic_resample( exp(weights), Np );
        pts = pts(parent);
        weights = log(ones(Np,1)/Np);
        postESS(kk) = Np;
        num_resam = num_resam + 1;
    else
        postESS(kk) = preESS(kk);
    end
    
    % Store last weights for next time
    last_weights = weights;
    
    % Store particles and weights
    for ii = 1:Np
        pts_store(ii).nonlin_samp(:,kk) = pts(ii).nonlin_samp(:,kk);
        pts_store(ii).lin_mn(:,kk) = pts(ii).lin_mn(:,kk);
        pts_store(ii).lin_vr(:,:,kk) = pts(ii).lin_vr(:,:,kk);
    end 
    wts_array{kk} = weights;
    
    % Store latest linear mean for filtering estimate
    est(kk) = sum(arrayfun(@(x,y) y*x.lin_mn(1,kk), pts, exp(last_weights)));

end

final_filt_pts = pts;

end

