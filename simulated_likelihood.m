function pk_hat = simulated_likelihood(observed_cells,k,dt,U_param,sigma,S,L)
% observed cells here is a cell array where cell i contains the current
% location of cell i
% propagate loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%observed_cells = calculate_neighbours(observed_cells,k);
N = length(observed_cells)-1;
pk_hat = cell(S,N,2);
dt = dt/L;
for s = 1:S
    copy = observed_cells;
    for l = 1:L-1
        F = interactions(observed_cells,U_param,k);
        for i = 1:N
        % quite slow but of ok structure
        
            if observed_cells{i}.d_time > k
                observed_cells{i}.location(1,k) = observed_cells{i}.location(1,k)...
                    + dt*F(1,i) + sqrt(dt)*sigma*randn(1,1);
                observed_cells{i}.location(2,k) = observed_cells{i}.location(2,k)...
                    + dt*F(2,i) + sqrt(dt)*sigma*randn(1,1);
                % addition of the randomness at this stage is ok since the
                % interactions have already been calculated.
            end
        end
    end
    F = interactions(observed_cells,U_param,k);
    for i = 1:N
        pk_hat{s,i,1} = dt*F(:,i) + observed_cells{i}.location(:,k); % the mean
        pk_hat{s,i,2} = eye(2)*sqrt(dt)*sigma;% the standard deviation
    end
    observed_cells = copy;
end
end