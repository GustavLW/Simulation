hold off
k1 = U_param(1);
a1 = U_param(2);
l1 = U_param(3);
k2 = U_param(4);
a2 = U_param(5);
l2 = U_param(6);

u   = @(r)  -a1 + a1./(1+exp(-k1*(r.^2-l1^2))) - a2*(1./(1+exp(-k2*(r.^2-l2^2))) - 1./(1+exp(-k1*(r.^2-l1^2))));
du  = @(r) 2*k1*(a1+a2)*r.*exp(k1*(r.^2+l1^2))./(exp(k1*l1^2)+exp(k1*r.^2)).^2 - ...
         2*k2*a2*r.*exp(k2*(r.^2+l2^2))./(exp(k2*l2^2)+exp(k2*r.^2)).^2;   
r = linspace(0,3,1001);

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
plot(r,-u(r))
hold on
grid on
%
plot(r,-du(r))
plot(r,-f(r))
%axis([0 3 -0.0025 0.0012])

%%

score = zeros(10,1);


R = 3;
t = 10;
for i = 1:10
    eps = 0.025*i;
    U_param = [7.0000    0.0016    0.3000    0.7500    0.0001    1.2500];
    sigma = 0.011;
    dt = observed_cells{end}*60;
    S   = 20;
    L   = 15;
    % rate agent is what we feed our agents into. contains evertthign!!!!
    ell = rate_agent(observed_cells,U_param,sigma,S,L);
    g = penalty_term(t,R,eps,U_param);
    score(i) = ell - g;
end
%%
%plot(0.0016*(0.5+(1:10)/10),score0)

r    = linspace(0,2,100001);
type = 1;
U_param = [0.0002 4];
v = u(r,type,U_param);
V = cumsum(v)*(r(2)-r(1));
V = V-V(end);
Q = v.*exp(-V);

subplot(3,1,1)
plot(r,v,'k')
grid on
axis([0 2 -0.2 0.2])
subplot(3,1,2)
plot(r,V,'k')
grid on
axis([0 2 -0.2 2])
subplot(3,1,3)

plot(r,Q,'k')
grid on
%%
plot(r,V,'k')
 grid on
axis([0 2 -0.001 0.005])


%%
a  = 4;
r0 = 1;
r = linspace(0,3,1001);
U = (1 - ((exp(r0)-1)./(exp(r) - 1))).^2 - 1;
plot(r,U,'k')
axis([0 3 -2 4])
grid on


