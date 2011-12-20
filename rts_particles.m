function [ pts ] = rts_particles( flags, params, pts, observs, init_mn, init_vr)
%KF_PARTICLES Run Kalman on particles

Np = length(pts);
P = params.ARO;
K = params.K;

if nargin > 3
    filter = true;
else
    filter = false;
end

for ii = 1:Np
    
    if filter
        forw_mn = zeros(P, K);
        forw_vr = zeros(P, P, K);
        forw_mn(:, 1) = init_mn;
        forw_vr(:, :, 1) = init_vr;
    else
        forw_mn = pts(ii).lin_mn;
        forw_vr = pts(ii).lin_vr;
    end
    
    A = zeros(P, P, K);
    Q = zeros(P, P, K);
    
    for kk = 2:K
        
        [A(:,:,kk-1), Q(:,:,kk-1)] = construct_transmats(flags, params, pts(ii).nonlin_samp(1:P,kk), pts(ii).nonlin_samp(P+1,kk));
        if filter
            H = [1 zeros(1, P-1)];
            R = params.noise_vr;
            [pred_mn, pred_vr] = kf_predict(forw_mn(:,kk-1), forw_vr(:,:,kk-1), A(:,:,kk-1), Q(:,:,kk-1));
            [forw_mn(:,kk), forw_vr(:,:,kk)] = kf_update_1D(pred_mn, pred_vr, observs(kk), H, R);
        end
        
    end
    
    back_mn = zeros(P,K);
    back_vr = zeros(P,P,K);
    [back_mn(:,P:end), back_vr(:,:,P:end)] = rts_smooth(forw_mn(:,P:end), forw_vr(:,:,P:end), A(:,:,P:end), Q(:,:,P:end));
    pts(ii).smooth_mn = back_mn;
    pts(ii).smooth_vr = back_vr;
    
end

end

