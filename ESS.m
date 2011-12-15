function [ ESS ] = ESS( w )
%ESS Calculate expected sample size

ESS = 1/( sum( exp(w).^2 ) );

end

