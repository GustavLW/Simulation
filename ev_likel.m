function ell = ev_likel(x,i,pk_hat)
ell = 0;
S = size(pk_hat,1);
for s = 1:S
    tmp_mu = pk_hat{s,i,1};
    tmp_sg = pk_hat{s,i,2};
    ell = ell + mvnpdf(x,tmp_mu',tmp_sg)/S;
end