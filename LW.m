%Solves the shallow water equations using the Lax-Wendroff scheme
function [H,m,E_h,E_m] = LW(dx, dt, T)
nx = 100/dx + 1;
Nt = T/dt;
[H0,m0,u0] = init(nx);
E1 = m0;
E2 = m0.*u0 + 0.5*9.81*H0.^2;
H(1:nx) = H0;
m(1:nx) = m0;
for j = 1:Nt
    H_i(1:nx-1) = 0.5*(H0(1:nx-1)+H0(2:nx)) - 0.5*dt/dx.*(E1(2:nx)-E1(1:nx-1));
    m_i(1:nx-1) = 0.5*(m0(1:nx-1)+m0(2:nx)) - 0.5*dt/dx.*(E2(2:nx)-E2(1:nx-1));
    E1_i = m_i;
    E2_i = m_i.^2./H_i + 0.5*9.81*H_i.^2;
    H(2:nx-1) = H0(2:nx-1) -dt/dx.*(E1_i(2:nx-1) - E1_i(1:nx-2));
    m(2:nx-1) = m0(2:nx-1) -dt/dx.*(E2_i(2:nx-1) - E2_i(1:nx-2));
    H0 = H;
    m0 = m;
    E1 = m;
    E2 = m.^2./H + 0.5*9.81*H.^2;
end
[H_a,m_a] = analytical(dx,T);
[E_h,E_m] = Error(H,H_a,m,m_a,nx);
end









