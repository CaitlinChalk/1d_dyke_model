%run model for simple example with unit pressure profile at the tip of the
%dyke

clear 
close all

%define spatial 1D grid
dz = 0.0236; 
z0 = -100;
z = (0:dz:11.8)';
n = length(z);


z2 = z + z0;
zf = z2(n) + dz/100;

%convert to frame of reference moving with the tip
x = zf - z2;

%Pressure profile that has an analytic solution
lambda = 2.36;
Pe0 = ones(n,1);
tf = (x > lambda);
Pe0(tf) = 0;

%% Calculate half width from pressure

%x = b/A => xA = b
%x = A\b => Ax = b
%h = M1*Pe => M1 = h/Pe, Pe = M1\h;
%Minv*h = Pe => Minv = Pe/h

%numerical solution
[M1,M2] = matrixM_new(x,dz,1,1);
h_numeric = M1*Pe0;

%analytic solution
h_analytic = analyticSolution(x,lambda);

close all
figure(1); hold on
plot(x,h_analytic,'-k')
plot(x,h_numeric,'o')
legend(["Analytic solution","Numerical solution"],"Interpreter","Latex")
xlabel("$z_f - z$","Interpreter","Latex")
ylabel("$h$","Interpreter","Latex")

%% Calculate pressure from half width

%recreate numeric pressure profile using \ operator
h2 = h_numeric;
Kc = 1;
h2(n) = Kc*(2^0.5);
Pe_numeric = M1\h2;

%recreate pressure profile using the inverse of M1
%estimate matrix inverse
M1inv = inv(M1); %using inv function
M1inv2 = real(Pe0)/real(h_numeric); %using / operator
Pe_numeric2 = M1inv*h2;

%recreate analytic pressure profile using iterative solver gmres
tol1 = 1e-06;
maxit = 10;
res = [];
Pe_gm = gmres(M1,h_analytic,res,tol1,maxit);

close all
figure(2); hold on
plot(x,Pe0,'-k') %analytical pressure profile
plot(x,Pe_numeric2,'o')
%plot(x,Pe_gm,'*')
legend(["Analytic solution","Numerical solution","gmres solution"],"Interpreter","Latex")
xlabel("$z_f - z$","Interpreter","Latex")
ylabel("$P_e$","Interpreter","Latex")



