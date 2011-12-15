function [A, Q] = construct_transmats(flags, params, ar_coefs, trans_noise)
%CONSTRUCT_TRANSMATS Construct the transition matrices for the linear state
%for progression using an AR model

P = length(ar_coefs);

A = [ar_coefs'; eye(P-1) zeros(P-1,1)];
Q = zeros(P); Q(1,1) = trans_noise;

% diag(fliplr(trans_noise));

end

