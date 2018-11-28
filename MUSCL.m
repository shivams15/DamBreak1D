%Solves the shallow water equations using Godunov scheme with MUSCL
%reconstruction of fluxes at the cell faces
function [H,m,E_h,E_m] = MUSCL(dx,dt,T)
nx = 100/dx + 1;
Nt = T/dt;
[H0,m0,u0] = init(nx);
H(1:nx) = H0;
m(1:nx) = m0;
S1(1:nx) = 0;
S2(1:nx) = 0;
E1_i(1:nx-1) = 0;
E2_i = E1_i;
for j = 1:Nt
    S1(1) = (H0(2) - H0(1))/dx;
    S2(1) = (m0(2) - m0(1))/dx;
    S1(nx) = (H0(nx) - H0(nx-1))/dx;
    S2(nx) = (m0(nx) - m0(nx-1))/dx;
    for i = 2:nx-1
        S1(i) = minmode((H0(i) - H0(i-1))/dx,(H0(i+1) - H0(i))/dx);
        S2(i) = minmode((m0(i) - m0(i-1))/dx,(m0(i+1) - m0(i))/dx);
    end
    for i = 1:nx-1
        H0_l = H0(i) + dx/2*S1(i);
        H0_r = H0(i+1) - dx/2*S1(i+1);
        m0_l = m0(i) + dx/2*S2(i);
        m0_r = m0(i+1) - dx/2*S2(i+1);
        alpha = abs(m0_l/H0_l + sqrt(9.81*H0_l));
        if alpha < abs(m0_r/H0_r + sqrt(9.81*H0_r))
            alpha = abs(m0_r/H0_r + sqrt(9.81*H0_r));
        end
        E1_i(i) = 0.5*(m0_l + m0_r) + 0.5*alpha*(H0_l - H0_r);
        E2_i(i) = 0.5*(m0_l^2/H0_l + 0.5*9.81*(H0_l^2 + H0_r^2) +m0_r^2/H0_r) + 0.5*alpha*(m0_l - m0_r);
    end
    H(2:nx-1) = H0(2:nx-1) - dt/dx.*(E1_i(2:nx-1) - E1_i(1:nx-2));
    m(2:nx-1) = m0(2:nx-1) - dt/dx.*(E2_i(2:nx-1) - E2_i(1:nx-2));
    H0 = H;
    m0 = m;
end
[H_a,m_a] = analytical(dx,T);
[E_h,E_m] = Error(H,H_a,m,m_a,nx);
end