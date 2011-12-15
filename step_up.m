function [ ar ] = step_up( ref )
%STEP_UP Levinson step-up recursion

[P, N] = size(ref);
ar = zeros(size(ref));

ar(1,:) = 1;

for jj = 0:P-1
    prev_ar = ar;
    for ii = 1:jj
        ar(ii,:) = prev_ar(ii,:) + ref(jj+1,:).*prev_ar(jj-ii+1,:);
    end
    ar(jj+1,:) = ref(jj+1,:);
end

ar = -ar;

end

