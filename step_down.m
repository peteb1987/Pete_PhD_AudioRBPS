function [ ref ] = step_down( ar )
%STEP_DOWN Levinson step-down recursion

ar = -ar;

[P, N] = size(ar);
ref = zeros(size(ar));

ref(P,:) = ar(P,:);

for jj = P-1:-1:1
    prev_ar = ar;
    for ii = 1:jj
        ar(ii,:) = (prev_ar(ii,:)-ref(jj+1,:).*prev_ar(jj-ii+1,:))./(1-ref(jj+1,:).^2);
    end
    ref(jj,:) = ar(jj,:);
end

end

