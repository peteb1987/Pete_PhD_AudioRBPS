function [ Nup, Nut ] = count_unique_particles( K, pts )
%COUNT_UNIQUE_PARTICLES Does what it says on the tin


Np = length(pts);
Nup = Np*ones(1,K);

% Loop through time
for k = 1:K
    
    if rem(k, 100)==0
        disp(['Time step ' num2str(k)]);
    end
    
    uniq = true(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        u_ii = pts(ii).nonlin_samp(:,k);
        
        % Loop through later particles
        for jj = ii+1:Np
            
            % Only look at particles which we still know to be unique
            if uniq(jj)
                
                u_jj = pts(jj).nonlin_samp(:,k);
                
                    % Compare jump sequences
                    if all(u_ii==u_jj)
                        
                        % Particle is not unique, cross it off
                        uniq(jj) = false;
                        Nup(k) = Nup(k) - 1;
                        
                    end
                
            end
            
        end
        
    end
    
end

Nut = Nup(end);

end

