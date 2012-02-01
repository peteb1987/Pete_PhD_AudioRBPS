function [ init_pts ] = initialise_particles( flags, params )
%INITIALISE_PARTICLES Set up particles for RB filter

% Sample reflection coefficients
init_ref = inf(params.ARO,params.Np);
for ii = 1:params.Np
    for pp = 1:params.ARO
        while abs(init_ref(pp,ii))>1
            init_ref(pp,ii) = normrnd(0, sqrt(params.init_ref_vr));
        end
    end
end

% Convert to AR coefficients and append a process variance term
init_nonlin = num2cell([ step_up(init_ref);
                     exp(normrnd(params.init_logprocvar_mn*ones(1,params.Np), sqrt(params.init_logprocvar_vr)))] , 1)';

% Create the particles
init_pts = struct('nonlin_samp', init_nonlin, 'lin_mn', zeros(params.ARO,1), 'lin_vr', 1E-20*eye(params.ARO));

end

