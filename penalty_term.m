function g = penalty_term(t,R,eps,U_param)

k1 = U_param(1);
a1 = U_param(2);
l1 = U_param(3);
k2 = U_param(4);
a2 = U_param(5);
l2 = U_param(6);


u   = @(r)  -a1 + a1./(1+exp(-k1*(r.^2-l1^2))) - a2*(1./(1+exp(-k2*(r.^2-l2^2))) - 1./(1+exp(-k1*(r.^2-l1^2))));
du  = @(r) 2*k1*(a1+a2)*r.*exp(k1*(r.^2+l1^2))./(exp(k1*l1^2)+exp(k1*r.^2)).^2 - ...
         2*k2*a2*r.*exp(k2*(r.^2+l2^2))./(exp(k2*l2^2)+exp(k2*r.^2)).^2;   

f1 = @(r) 4*a1*k1^2*r.^2.*exp(k1*(l1^2+r.^2))./(exp(k1*l1^2+exp(k1*r.^2))).^2;
f2 = @(r) 8*a1*k1^2*r.^2.*exp(k1*(l1^2+r.^2)+k1*r.^2)./(exp(k1*l1^2+exp(k1*r.^2))).^3;
f3 = @(r) 2*a1*k1*exp(k1*(l1^2+r.^2))./(exp(k1*l1^2+exp(k1*r.^2))).^2;
f4 = @(r) 4*a2*k1^2*r.^2.*exp(k1*(l1^2+r.^2))./(exp(k1*l1^2+exp(k1*r.^2))).^2;
f5 = @(r) 8*a2*k1^2*r.^2.*exp(k1*(l1^2+r.^2)+k1*r.^2)./(exp(k1*l1^2+exp(k1*r.^2))).^3;
f6 = @(r) 2*a2*k1*exp(k1*(l1^2+r.^2))./(exp(k1*l1^2+exp(k1*r.^2))).^2;
f7 = @(r) 4*a2*k2^2*r.^2.*exp(k2*(l2^2+r.^2))./(exp(k2*l2^2+exp(k2*r.^2))).^2;
f8 = @(r) 8*a2*k2^2*r.^2.*exp(k2*(l2^2+r.^2)+k2*r.^2)./(exp(k2*l2^2+exp(k2*r.^2))).^3;
f9 = @(r) 2*a2*k2*exp(k2*(l2^2+r.^2))./(exp(k2*l2^2+exp(k2*r.^2))).^2;

f = @(r) f1(r)-f2(r)+f3(r)+f4(r)-f5(r)+f6(r)-f7(r)+f8(r)-f9(r);





g1 = max(0,abs(u(R)) - eps*abs(u(1)))^2
g2 = max(0,-u(0))^2;
g3 = max(0,du(0))^2;
g4 = max(0,u(1))^2;
g5 = max(0,du(1))^2;
g6 = max(0,-du(1))^2;
g7 = max(0,-f(1))^2;

g = 2^t*(g1 + g2 + g3 + g4 + g5 + g6 + g7);