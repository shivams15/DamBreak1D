function [H,m,u] = init(nx)
H(1:(nx+1)/2) = 10;
H((nx+3)/2:nx) = 1;
m(1:nx) = 0;
u(1:nx) = 0;
end