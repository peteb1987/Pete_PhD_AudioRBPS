function [ est ] = rts_particles( flags, params, pts, weights, observs )
%KF_PARTICLES Run Kalman on particles

Np = length(pts);
P = params.ARO;
K = params.K;

pts_smooth_mn = zeros(P,K,Np);

for ii = 1:Np
    
    forw_mn = zeros(P, K);
    forw_vr = zeros(P, P, K);
    
    mn = zeros(P,1);
    vr = zeros(P,P);
    
    A = zeros(P, P, K);
    Q = zeros(P, P, K);
    
    for kk = 1:K
        
        % Kalman filter
        [A(:,:,kk), Q(:,:,kk)] = construct_transmats(flags, params, pts(ii).nonlin_samp(1:P,kk), pts(ii).nonlin_samp(P+1,kk));
        H = [1 zeros(1, P-1)];
        R = params.noise_vr;
        [mn, vr] = kf_predict(mn, vr, A(:,:,kk), Q(:,:,kk));
        [mn, vr] = kf_update_1D(mn, vr, observs(kk), H, R);
        
        % Store
        forw_mn(:,kk) = mn;
        forw_vr(:,:,kk) = vr;
        
    end
    
    back_mn = zeros(P,K); back_vr = zeros(P,P,K);
    [back_mn(:,P:end), back_vr(:,:,P:end)] = rts_smooth(forw_mn(:,P:end), forw_vr(:,:,P:end), cat(3,A(:,:,P+1:end),zeros(P)), cat(3,Q(:,:,P+1:end),zeros(P)));
%     pts(ii).smooth_mn = back_mn;
%     pts(ii).smooth_vr = back_vr;
    
    pts_smooth_mn(:,:,ii) = back_mn;

end

est = sum( bsxfun(@times, permute(pts_smooth_mn,[3,1,2]), exp(weights)) , 1);
est = squeeze(est(1, 1, :));

end

