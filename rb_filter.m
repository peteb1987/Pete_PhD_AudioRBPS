function [ lin_mn_est, pts ] = rb_filter( flags, params, init_pts, observs )
%RB_FILTER Rao-Blackwellised particle filter

Np = params.Np;

% Initialise some things
pts = init_pts;
K = length(observs);
weights = cell(K,1);
last_weights = zeros(Np, 1);
lin_mn_est = zeros(size(observs));

for ii = 1:Np
    pts(ii).lin_vr = cat(3, pts(ii).lin_vr, repmat(1E-20*eye(params.ARO), [1 1 K-1]));
    pts(ii).lin_mn = [pts(ii).lin_mn, repmat(zeros(params.ARO,1), [1 K-1])];
    pts(ii).nonlin_samp = [pts(ii).nonlin_samp, repmat(zeros(params.ARO+1,1), [1 K-1])];
end
    
tic;

% Loop through time
for kk = 1:K
    
    if mod(kk,100)==0
        fprintf(1, '*** Time point %u. Last 100 took %f s.\n', kk, toc)
        tic;
    end
    
    weights{kk} = zeros(Np, 1);
    
    % Get the index of the last time. This ensures that the time 0 prior is overwritten
    last_idx = max(1,kk-1);
    P = min(params.ARO, kk);
    
    % Loop through particles
    for ii = 1:Np
        
        prev_nonlin_samp = pts(ii).nonlin_samp(:,last_idx);
        prev_lin_mn = pts(ii).lin_mn(1:P,last_idx);
        prev_lin_vr = pts(ii).lin_vr(1:P,1:P,last_idx);
%         pts(ii).lin_vr(:,:,kk) = 1E-20*eye(params.ARO);
        
        
        % Propose new nonlinear state (from transition density)
        pts(ii).nonlin_samp(:,kk) = sample_nonlin_transdens(flags, params, prev_nonlin_samp);
        
        % Run Kalman filter
        [A, Q] = construct_transmats(flags, params, pts(ii).nonlin_samp(1:P,kk), pts(ii).nonlin_samp(params.ARO+1,kk));
        [pred_mn, pred_vr] = kf_predict(prev_lin_mn, prev_lin_vr, A, Q);
        assert(isposdef(pred_vr));
        H = [1 zeros(1, P-1)];
        R = params.noise_vr;
        [pts(ii).lin_mn(1:P,kk), pts(ii).lin_vr(1:P,1:P,kk), ~, ~, ~, pred_lhood] = kf_update_1D(pred_mn, pred_vr, observs(kk), H, R);
%         H = eye(P); R = params.noise_vr*eye(P);
%         [pts(ii).lin_mn(1:P,kk), pts(ii).lin_vr(1:P,1:P,kk), ~, ~, ~, pred_lhood] = kf_update(pred_mn, pred_vr, observs(kk-P+1:kk), H, R);
        pts(ii).lin_vr(1:P,1:P,kk) = (pts(ii).lin_vr(1:P,1:P,kk)+pts(ii).lin_vr(1:P,1:P,kk)')/2;
        assert(isreal(pred_lhood));
        assert(isposdef(pts(ii).lin_vr(1:P,1:P,kk)));
        
        % Update weight
        weights{kk}(ii) = last_weights(ii) + log(pred_lhood);
        
    end
    
    % Normalise weights
    lin_weights = exp(weights{kk}); lin_weights = lin_weights/sum(lin_weights);
    weights{kk} = log(lin_weights);
    
    
    % Systematic resampling
    if ESS(weights{kk})<0.5*Np
        [ ~, parent ] = systematic_resample( exp(weights{kk}), Np );
        pts = pts(parent);
        weights{kk} = log(ones(Np,1)/Np);
    end
    
    % Store last weights for next time
    last_weights = weights{kk};
    
%     lin_mn_est = sum(multiprod(cat(3, pts.lin_mn), exp(last_weights), 3, 1), 3);
    lin_mn_est(kk) = sum(arrayfun(@(x,y) y*x.lin_mn(1,kk), pts, exp(last_weights)));
    
end

% lin_mn_est = [lin_mn_est(params.ARO+1:end); zeros(params.ARO,1)];

end

