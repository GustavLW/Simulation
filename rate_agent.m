function ell = rate_agent(observed_cells,U_param,sigma,S,L)

K   = size(observed_cells{1}.location,2);
dt = observed_cells{end}*60;
ell = 0;

for k = 1:K-1
    pk_hat = simulated_likelihood(observed_cells,k,dt,U_param,sigma,S,L);
    for i = 1:size(observed_cells,1)-1
        pik    = log(ev_likel(observed_cells{i}.location(:,k+1)',i,pk_hat));
        ell = ell + pik;
    end
end
