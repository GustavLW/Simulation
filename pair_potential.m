function f = pair_potential(xi,xj,U_param) 
% calculate the influence of cell j on cell i 
type = 1; % 0 for LJ, 1 for Morse, 2 for special sauce 
r    = norm(xi-xj);
if r < 3 % extremely arbitrary
    f = -(xi-xj)*u(r,type,U_param)/r;
else
    f = [0;0];
end
