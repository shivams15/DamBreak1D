function PlotSolutionMUSCL(dx,dt,T)
[H1,m1,~,~] = MUSCL(dx,dt,T);
figure(1)
plot(0:dx:100,H1);
xlabel('x (m)');
ylabel('H (m)');

figure(2)
plot(0:dx:100,m1./H1);
xlabel('x (m)');
ylabel('U (m/s)');

end