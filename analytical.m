%Generates the analytical solution for the 1D dam break problem
function [H,m] = analytical(dx,T)
nx = 100/dx + 1;
us = 7.35196;
hs = 3.96173;
s = 9.8193;
g = 9.81;
c1 = sqrt(g*10);
for i = 1:nx
    x = (i-1)*dx - 50;
    if x <= - c1*T
        H(i) = 10;
        u(i) = 0;
    elseif x <= T*(1.5*us - c1)
        H(i) = (2*c1 - x/T)^2/(9*g);
        u(i) = 2*(c1 + x/T)/3;
    elseif x <= s*T
        H(i) = hs;
        u(i) = us;
    else
        H(i) = 1;
        u(i) = 0;
    end
end
m = H.*u;
end