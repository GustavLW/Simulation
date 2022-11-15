function v = u(r,type,U_param)
% this function simply spits out the signed magnitude of the potential at a
% certain distance, without taking orientation into account. type = 0 for
% LJ and 1 for Morse; i will probably add more later berhaps.


if type == 0
    Vmin = U_param(1);
    alpha = U_param(2);
    D     = -Vmin/(2^(-2)-2^(-1));
    rho   = 2^(-1/alpha);
    v = -(D*alpha./r)*(2*(rho./r).^(2*alpha)-(rho./r).^(alpha));
elseif type == 1
    Vmin = U_param(1);
    alpha = U_param(2);
    v = -2*Vmin*alpha*(exp(-2*alpha*(r-1))-exp(-alpha*(r-1)));
elseif type == 2
    k1 = U_param(1);
    a1 = U_param(2);
    l1 = U_param(3);
    k2 = U_param(4);
    a2 = U_param(5);
    l2 = U_param(6);
    v  = 2*k1*r*(a1+a2)*exp(k1*(r^2+l1^2))/(exp(k1*l1^2)+exp(k1*r^2))^2 - ...
         2*k2*r*a2*exp(k2*(r^2+l2^2))/(exp(k2*l2^2)+exp(k2*r^2))^2;
else
    disp('Mata in en vettig potential!')
    v = [];
end