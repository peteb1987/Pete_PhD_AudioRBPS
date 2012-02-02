function [ data, ar, proc_var ] = generate_data( flags, params )
%GENERATE_DATA

K = params.K;
P = params.ARO;

data = zeros(K,1);
ar = zeros(P,K);
proc_var = zeros(1,K);

% Initialse ar coefficients
ref = inf(P,1);
for pp = 1:P
    while abs(ref(pp)) > 1
        ref(pp) = mvnrnd(0, params.init_ref_vr);
    end
end
ar(:,1) = step_up(ref);

% Initialise noise_var
log_proc_var = -inf;
while log_proc_var < params.min_log_proc_var
    log_proc_var = mvnrnd(params.init_logprocvar_mn, params.init_logprocvar_vr);
end
proc_var(1) = exp(log_proc_var);

% Loop through time
for kk = 2:K
    
    P = min(kk-1, params.ARO);
    
    % Propagate ar coefficients
    last_ref = step_down(ar(:,kk-1));
    ref = inf(size(last_ref));
    for pp = 1:params.ARO
        while abs(ref(pp)) > 1
            ref(pp) = mvnrnd(last_ref(pp), params.ref_trans_vr);
        end
    end
    ar(:,kk) = step_up(ref);
    
    % Propagate noise variance
    log_proc_var_mn = log(params.procvar_decay*proc_var(kk-1));
    log_proc_var = -inf;
    while log_proc_var < params.min_log_proc_var
        log_proc_var = mvnrnd(log_proc_var_mn, params.logprocvar_vr);
    end
    proc_var(kk) = exp(log_proc_var);
    
    % Sample linear state
    lin_mn = ar(1:P,kk)'*data(kk-1:-1:kk-P);
    data(kk) = mvnrnd(lin_mn, proc_var(kk));
    
end
    
end
