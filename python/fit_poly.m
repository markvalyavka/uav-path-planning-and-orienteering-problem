clear all
close all
samples = table2array(readtable('../samples_equidistant_08.csv','NumHeaderLines',1)); 
poly = table2array(readtable('../polynomials.csv','NumHeaderLines',1)); 
tsam = samples(:,1);
xsam = samples(:,2);
ysam = samples(:,3);
zsam = samples(:,4);

xvsam = samples(:,9);
yvsam = samples(:,10);
zvsam = samples(:,11);

xasam = samples(:,15);
yasam = samples(:,16);
zasam = samples(:,17);

plot3(xsam,ysam,zsam,'*')
hold on

p = quadprog()

P = [];
V = [];
A = [];
J = [];
NA = [];
T = [];
for i=1:3:size(poly)
   xpoly =  poly(i,3:end);
   ypoly =  poly(i+1,3:end);
   zpoly =  poly(i+2,3:end);
   tfrom =xpoly(1)
   tto =xpoly(2)-tfrom
   for t=0:0.001:tto 
       tau = [t^5,t^4,t^3,t^2,t,1];
       tauv = [5*t^4,4*t^3,3*t^2,2*t,1,0];
       taua = [20*t^3,12*t^2,6*t,2,0,0];
       tauj = [60*t^2,24*t,6,0,0,0];
       
       x = xpoly(3:end);
       xp = tau*x';
       xv = tauv*x';
       xa = taua*x';
       xj = tauj*x';
       y = ypoly(3:end);
       yp = tau*y';
       yv = tauv*y';
       ya = taua*y';
       yj = tauj*y';
       z = zpoly(3:end);
       zp = tau*z';
       zv = tauv*z';
       za = taua*z';
       zj = tauj*z';
       P = [P;[xp,yp,zp]];
       V = [V;[xv,yv,zv]];
       A = [A;[xa,ya,za]];
       J = [J;[xj,yj,zj]];
       NA = [NA;norm([xa,ya,za])];
       T = [T;t+tfrom];
   end
end




plot3(P(:,1),P(:,2),P(:,3))


figure(2)
hold on
plot(T,V(:,1))
plot(tsam,xvsam,'.--')
plot(T,V(:,2))
plot(tsam,yvsam,'.--')
plot(T,V(:,3))
plot(tsam,zvsam,'.--')


figure(3)
hold on
plot(T,A(:,1))
plot(tsam,xasam,'.--')
plot(T,A(:,2))
plot(tsam,yasam,'.--')
plot(T,A(:,3))
plot(tsam,zasam,'.--')
plot(T,NA,'c--')

figure(4)
hold on
plot(T,J(:,1))
plot(T,J(:,2))
plot(T,J(:,3))
ylim([-320,320])

