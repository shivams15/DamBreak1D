%Calculates the error in the computed solution
function [E_h,E_m] = Error(H,H_a,m,m_a,N)
diff = H - H_a;
E_h = sqrt(sum(diff.*diff))/N;
diff = m - m_a;
E_m = sqrt(sum(diff.*diff))/N;
end