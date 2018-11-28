%Solves the shallow water equations using the Godunov scheme with
%Lax-Friedrich flux
function [H,m,E_h,E_m] = Godunov(dx,dt,T)
nx = 100/dx + 1;
Nt = T/dt;
[H0,m0,u0] = init(nx);
E1 = m0;
E2 = m0.*u0 + 0.5*9.81*H0.^2;
H(1:nx) = H0;
m(1:nx) = m0;
E1_i(1:nx-1) = 0;
E2_i = E1_i;
for j = 1:Nt
    for i = 1:nx-1
        alpha = abs(u0(i) + sqrt(9.81*H0(i)));
        if alpha < abs(u0(i+1) + sqrt(9.81*H0(i+1)))
            alpha = abs(u0(i+1) + sqrt(9.81*H0(i+1)));
        end
        E1_i(i) = 0.5*(E1(i) + E1(i+1)) + 0.5*alpha*(H0(i) - H0(i+1));
        E2_i(i) = 0.5*(E2(i) + E2(i+1)) + 0.5*alpha*(m0(i) - m0(i+1));
    end
    H(2:nx-1) = H0(2:nx-1) - dt/dx.*(E1_i(2:nx-1) - E1_i(1:nx-2));
    m(2:nx-1) = m0(2:nx-1) - dt/dx.*(E2_i(2:nx-1) - E2_i(1:nx-2));
    H0 = H;
    m0 = m;
    u0 = m0./H0;
    E1 = m;
    E2 = m.^2./H + 0.5*9.81*H.^2;
end
[H_a,m_a] = analytical(dx,T);
[E_h,E_m] = Error(H,H_a,m,m_a,nx);
end