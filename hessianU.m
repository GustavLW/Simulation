function hU = hessianU(xik,xjk,type,U_param)
if type == 2
    k1 = U_param(1);
    a1 = U_param(2);
    l1 = U_param(3);
    k2 = U_param(4);
    a2 = U_param(5);
    l2 = U_param(6);
    
    r = xik-xjk;
    R = norm(r);
    
    g1 = -k1*(R^2-l1^2);
    g2 = -k4*(R^2-l2^2);
    L1 = exp(g1)/(1+exp(g1))^2;
    L2 = exp(g2)/(1+exp(g2))^2;
    G11 = -2*k1*r(1);
    G12 = -2*k1*r(2);
    G21 = -2*k2*r(1);
    G22 = -2*k2*r(2);
    A1  = exp(2*g1)/(1+exp(g1))^3;
    A2  = exp(2*g2)/(1+exp(g2))^3;
    
    p11 = -L1*G11^2 + 2*A1*G11+2*k1*L1;
    p12 = -L1*G11*G12 + 2*A1*G11*G12;
    p22 = -L1*G12^2 + 2*A1*G12 + 2*k1*L1;
    
    q11 = -L2*G21^2 + 2*A2*G21+2*k2*L2;
    q12 = -L2*G21*G22 + 2*A2*G21*G22;
    q22 = -L2*G22^2 + 2*A2*G22 + 2*k2*L2;
    
    hU  = [a1*p11-a2*(q11-p11), a1*p12-a2*(q12-p12);...
           a1*p12-a2*(q12-p12), a1*p22-a2*(q22-p22)];
end