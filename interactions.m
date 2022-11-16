function F = interactions(population,U_param,k)
if size(population{1}.neighbours,1) == 1
    N_nei = size(population{1}.neighbours,2);
    N     = length(population);
    x     = zeros(2,N);
    F     = zeros(2,N);
    for i = 1:N
        x(:,i) = population{i}.location;
    end
    
    for i = 1:N
        xi  = x(:,i);
        nei = population{i}.neighbours;
        if sum(size(nei)) > 0
            for j = 1:min(N,N_nei)
                if nei(j) == 1
                    xj     = x(:,j);
                    f      = pair_potential(xi,xj,U_param);
                    F(:,i) = F(:,i) + f;
                end
            end
        end
    end
else
    N_nei = size(population{1}.neighbours,2);
    N     = length(population)-1;
    x     = zeros(2,N);
    F     = zeros(2,N);
    for i = 1:N
        x(:,i) = population{i}.location(:,k);
    end
    
    for i = 1:N
        xi  = x(:,i);
        nei = population{i}.neighbours(k,:);
        if sum(size(nei)) > 0
            for j = 1:min(N,N_nei)
                if nei(j) == 1
                    xj = x(:,j);
                    F(:,i) = F(:,i) + pair_potential(xi,xj,U_param);
                end
            end
        end
    end
end
