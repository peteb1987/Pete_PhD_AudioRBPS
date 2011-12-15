function [ ar ] = step_up_check( ref )
%STEP_UP Levinson step-up recursion

[P, ~] = size(ref);
ar = zeros(P+1,P+1);

ar(1,1) = 1;

for nn = 1:P
    n_ind = nn+1;
    for ii = 0:nn-1
        i_ind = ii+1;
        ar(n_ind,i_ind) = ar(n_ind-1,i_ind) + ref(nn)*ar(n_ind-1,nn-ii+1);
    end
    ar(n_ind,n_ind) = ref(nn);
end


end

