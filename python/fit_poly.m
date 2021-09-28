clear all
close all
samples = table2array(readtable('../samples_equidistant_08.csv','NumHeaderLines',1)); 
poly = table2array(readtable('../polynomials.csv','NumHeaderLines',1)); 
t = samples(:,1);
x = samples(:,2);
y = samples(:,3);
z = samples(:,4);

plot3(x,y,z,'*')
hold on

P = []
for i=1:3:size(poly)
   xpoly =  poly(i,3:end);
   ypoly =  poly(i+1,3:end);
   zpoly =  poly(i+2,3:end);
   tfrom =xpoly(1)
   tto =xpoly(2)-tfrom
   for t=0:0.01:tto 
       tau = [t^5,t^4,t^3,t^2,t,1];
       
       x = xpoly(3:end)
       xp = tau*x';
       y = ypoly(3:end);
       yp = tau*y';
       z = zpoly(3:end);
       zp = tau*z';
       P = [P;[xp,yp,zp]]    ;
   end
   
   
end


plot3(P(:,1),P(:,2),P(:,3))


