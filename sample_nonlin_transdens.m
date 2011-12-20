function [ new_nonlin_samp, prob ] = sample_nonlin_transdens( flags, params, old_nonlin_samp, new_nonlin_samp )
%SAMPLE_NONLIN_TRANSDENS Sample transition density for audio model

old_ar = old_nonlin_samp(1:end-1, 1);
old_proc_vr = old_nonlin_samp(end, 1);

% Convert to reflection coefficients
old_ref = step_down(old_ar);

% Log process variance mean
log_change_mn = log(params.procvar_decay*old_proc_vr);

if nargin < 4
    
    % Rejection sample from truncated Gaussian
    new_ref = inf(size(old_ref));
    for jj = 1:length(new_ref)
        while abs(new_ref(jj))>1
            new_ref(jj) = normrnd(old_ref(jj), sqrt(params.ref_trans_vr));
        end
    end
    
    % Convert back to ar coefficients
    new_ar = step_up(new_ref);
    
    % Sample process var
    log_change = -inf;
    while log_change_mn + log_change < -10
        log_change = normrnd(0, sqrt(params.logprocvar_vr));
        new_proc_vr = exp( log_change_mn + log_change );
    end
    assert(isreal(new_proc_vr));
    
%     % LOWER LIMIT THE PROCESS VARIANCE - FUDGE
%     new_proc_vr = max(new_proc_vr, 1E-8);
%     if new_proc_vr<1E-15
%         warning('Teeny Weeny Process Variance!!');
%     end
    
    % Concatenate nonlinear state
    new_nonlin_samp = [new_ar; new_proc_vr];
    
else
    
    new_ar = new_nonlin_samp(1:end-1, 1);
    new_proc_vr = new_nonlin_samp(end, 1);
    log_change = log(new_proc_vr/(old_proc_vr*params.procvar_decay));
    
    % Convert to reflection coefficients
    new_ref = step_down(new_ar);
    
end

if nargout > 1
    
    % Calculate probability
    procvar_prob = log(normpdf(log_change, 0, sqrt(params.logprocvar_vr)));
%     for jj = 1:length(new_ref)
%         prob = prob + log( normpdf(new_ref(jj), old_ref(jj), sqrt(params.ref_trans_vr)) );
%         Z = normcdf(1, old_ref(jj), sqrt(params.ref_trans_vr)) - normcdf(-1, old_ref(jj), sqrt(params.ref_trans_vr));
%         prob = prob - log(Z);
%     end
    ref_prob = log( normpdf(new_ref, old_ref, sqrt(params.ref_trans_vr)) );
    norm_consts = log( normcdf(ones(size(old_ref)), old_ref, sqrt(params.ref_trans_vr)) - normcdf(-ones(size(old_ref)), old_ref, sqrt(params.ref_trans_vr)) );
    prob = procvar_prob + sum(ref_prob) - sum(norm_consts);
    
end

if any(abs(roots([1, -new_ar']))>1)
    warning('Not stable!');
end

end

