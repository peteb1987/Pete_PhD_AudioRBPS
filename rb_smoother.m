function [ smooth_pts ] = rb_smoother( flags, params, filt_pts, filt_weights, observs )
%RB_SMOOTHER Run Rao-Balckwellised particle smoother

Ns = params.Ns;
Np = params.Np;
K = length(observs);
lin_mn_est = zeros(size(observs));

smooth_pts = filt_pts(1:Ns);

% Loop through smoothing trajectories
for ss = 1:Ns
    
    ss
    
    smooth_pts(ss).lin_mn = zeros(params.ARO, K);
    smooth_pts(ss).lin_vr = zeros(params.ARO, params.ARO, K);
    smooth_pts(ss).lin_mn(:,K) = smooth_pts(ss).lin_mn(:,K);
    smooth_pts(ss).lin_vr(:,:,K) = 1E-3*eye(params.ARO);
    
    % Loop backwards in time
    for kk = K-params.ARO:-1:params.ARO
        
        if mod(kk,100)==0
            fprintf(1, '*** Time point %u.\n', kk)
        end
        
        P = min(params.ARO, kk);
        
        % Calculate the sampling weights
        sampling_weights = zeros(Np, 1);
        trans_prob = zeros(Np, 1);
        lin_prob = zeros(Np, 1);
        
        % Loop through the previous frame particles
        for jj = 1:Np
            
            % Calculate nonlinear transition density
            [~, trans_prob(jj)] = sample_nonlin_transdens(flags, params, filt_pts(jj).nonlin_samp(:,kk), smooth_pts(ss).nonlin_samp(:,kk+1));
%             sampling_weights(jj) = sampling_weights(jj) + trans_prob;
            
%             % Kalman smoother recursions - P-step backward prediction
%             pred_mn = smooth_pts(ss).lin_mn(1:P,kk+P);
%             pred_vr = smooth_pts(ss).lin_vr(1:P,1:P,kk+P);
%             sum_log_det_A = 0;
%             for pp = P:-1:1
%                 [A, Q] = construct_transmats(flags, params, filt_pts(jj).nonlin_samp(1:P,kk+pp), filt_pts(jj).nonlin_samp(params.ARO+1,kk+pp));
%                 [pred_mn, pred_vr] = kf_predict(pred_mn, pred_vr, inv(A), (A\Q/A'));
%                 sum_log_det_A = sum_log_det_A + log(abs(det(A)));
%             end
%             pred_vr = (pred_vr+pred_vr')/2;
% %             assert(isposdef(pred_vr));
% %             assert(isposdef(filt_pts(jj).lin_vr(1:P,1:P,kk)+pred_vr));
%             
%             test_vr = filt_pts(jj).lin_vr(1:P,1:P,kk)+pred_vr;
%             test_mn = filt_pts(jj).lin_mn(1:P,kk) - pred_mn;
%             lin_prob(jj) = log( mvnpdf( zeros(1,P), test_mn', test_vr) ) + sum_log_det_A;

            % P-step forward prediction
            pred_mn = filt_pts(jj).lin_mn(1:P,kk);
            pred_vr = filt_pts(ss).lin_vr(1:P,1:P,kk);
            for pp = 1:P
                [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+pp), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+pp));
                [pred_mn, pred_vr] = kf_predict(pred_mn, pred_vr, A, Q);
            end
            pred_vr = (pred_vr+pred_vr')/2;
            
            % Linear probability
            test_vr = smooth_pts(ss).lin_vr(1:P,1:P,kk+P)+pred_vr;
            test_mn = smooth_pts(ss).lin_mn(1:P,kk+P) - pred_mn;
            lin_prob(jj) = log( mvnpdf( zeros(1,P), test_mn', test_vr) );
            
            % Add up weights
            sampling_weights(jj) = filt_weights{kk}(jj) ...
                + trans_prob(jj) ...
                + lin_prob(jj);
            
        end
        
        % Backwards sampling
        samp_ind = randsample(Np, 1, true, exp(sampling_weights));
        smooth_pts(ss).nonlin_samp(:, 1:kk) = filt_pts(samp_ind).nonlin_samp(:, 1:kk);
        
        % Backwards Kalman smoother update
        [A, Q] = construct_transmats(flags, params, smooth_pts(ss).nonlin_samp(1:P,kk+P), smooth_pts(ss).nonlin_samp(params.ARO+1,kk+P));
        [pred_mn, pred_vr] = kf_predict(smooth_pts(ss).lin_mn(1:P,kk+P), smooth_pts(ss).lin_vr(1:P,1:P,kk+P), inv(A), A\Q/A');
        H = [zeros(1, P-1) 1];
        R = params.noise_vr;
        [smooth_pts(ss).lin_mn(1:P,kk+P-1), smooth_pts(ss).lin_vr(1:P,1:P,kk+P-1)] = kf_update_1D(pred_mn, pred_vr, observs(kk), H, R);
%         H = eye(P); R = params.noise_vr*eye(P);
%         [smooth_pts(ss).lin_mn(1:P,kk-1), smooth_pts(ss).lin_vr(1:P,1:P,kk-1)] = kf_update(pred_mn, pred_vr, observs(kk-P:kk-1), H, R);
        smooth_pts(ss).lin_vr(1:P,1:P,kk+P-1) = (smooth_pts(ss).lin_vr(1:P,1:P,kk+P-1)+smooth_pts(ss).lin_vr(1:P,1:P,kk+P-1)')/2;
        assert(isposdef(smooth_pts(ss).lin_vr(1:P,1:P,kk+P-1)));
        
    end
    
end




end

