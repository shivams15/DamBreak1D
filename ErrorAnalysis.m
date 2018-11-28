figure(1)
dx = [0.05 0.125 0.25];
for z = 1:3
    [~,~,e1_h(z),e1_m(z)] = Godunov(dx(z),0.0005,4);
    [~,~,e2_h(z),e2_m(z)] = MUSCL(dx(z),0.0005,4);
    [~,~,e3_h(z),e3_m(z)] = LW(dx(z),0.0005,4);
end
plot(log(dx),log(e1_m),log(dx),log(e2_m),log(dx),log(e3_m));
xlabel('log({\Delta}x)');
ylabel('log(error)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff');

figure(2)
plot(log(dx),log(e1_h),log(dx),log(e2_h),log(dx),log(e3_h));
xlabel('log({\Delta}x)');
ylabel('log(error)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff');


figure(3)
dt = [0.0005 0.001 0.002];
for z = 1:3
    [~,~,e1_h(z),e1_m(z)] = Godunov(0.05,dt(z),4);
    [~,~,e2_h(z),e2_m(z)] = MUSCL(0.05,dt(z),4);
    [~,~,e3_h(z),e3_m(z)] = LW(0.05,dt(z),4);
end
plot(log(dt),log(e1_m),log(dt),log(e2_m),log(dt),log(e3_m));
xlabel('log({\Delta}t)');
ylabel('log(error)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff');

figure(4)
plot(log(dt),log(e1_h),log(dt),log(e2_h),log(dt),log(e3_h));
xlabel('log({\Delta}t)');
ylabel('log(error)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff');

[H_a,m_a] = analytical(0.05,4);

figure(5)
[H1,m1,~,~] = Godunov(0.05,0.0005,4);
plot(0:0.05:100,H1,0:0.05:100,H_a,':');
xlabel('x (m)');
ylabel('H (m)');
legend('Godunov 1st order', 'Analytical');

figure(6)
plot(0:0.05:100,m1./H1,0:0.05:100,m_a./H_a,':');
xlabel('x (m)');
ylabel('U (m/s)');
legend('Godunov 1st order', 'Analytical');

figure(7)
[H2,m2,~,~] = MUSCL(0.05,0.0005,4);
plot(0:0.05:100,H2,0:0.05:100,H_a,':');
xlabel('x (m)');
ylabel('H (m)');
legend('Godunov with MUSCL', 'Analytical');

figure(8)
plot(0:0.05:100,m2./H2,0:0.05:100,m_a./H_a,':');
xlabel('x (m)');
ylabel('U (m/s)');
legend('Godunov with MUSCL', 'Analytical');

figure(9)
[H3,m3,~,~] = LW(0.05,0.002,4);
plot(0:0.05:100,H3,0:0.05:100,H_a,':');
xlabel('x (m)');
ylabel('H (m)');
legend('Lax-Wendroff', 'Analytical');

figure(10)
plot(0:0.05:100,m3./H3,0:0.05:100,m_a./H_a,':');
xlabel('x (m)');
ylabel('U (m/s)');
legend('Lax-Wendroff', 'Analytical');

figure(17)
plot(0:0.05:100,H1,0:0.05:100,H2,0:0.05:100,H3,0:0.05:100,H_a);
xlabel('x (m)');
ylabel('H (m)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff', 'Analytical');

figure(18)
plot(0:0.05:100,m1./H1,0:0.05:100,m2./H2,0:0.05:100,m3./H3,0:0.05:100,m_a./H_a);
xlabel('x (m)');
ylabel('U (m/s)');
legend('Godunov 1st order','Godunov with MUSCL','Lax-Wendroff', 'Analytical');

figure(11)
dt = [0.001 0.0025];
for z = 1:2
    [H,~,~,~] = MUSCL(0.05,dt(z),4);
    plot(0:0.05:100,H);
    hold on;
end 

figure(12)
dt = [0.001 0.0025];
for z = 1:2
    [H,m,~,~] = MUSCL(0.05,dt(z),4);
    plot(0:0.05:100,m./H);
    hold on;
end 

figure(13)
dt = [0.0005 0.0025];
for z = 1:2
    [H,~,~,~] = LW(0.05,dt(z),4);
    plot(0:0.05:100,H);
    hold on;
end 

figure(14)
dt = [0.0005 0.0025];
for z = 1:2
    [H,m,~,~] = LW(0.05,dt(z),4);
    plot(0:0.05:100,m./H);
    hold on;
end 



