%script to call dykeModel functions and solve system
clear 
close all

mu = 0.01; %viscosity (Pa s)
nu = 0.2;  %Poisson ratio
G = 10^10; %Shear modulus (Pa)
rhom = 2700; %Magma density (kg /m^3) 
rhos = 2900; %Host rock density
drho = rhos - rhom;
g = 9.81;
Q = 1;

%define input parameters to pass to model
parameters = struct();

%numerical parameters
parameters.dt = 0.000001;
parameters.dz = 0.0236; 
parameters.nTimeSteps = 5000000;
parameters.tol = 1.e-7;
parameters.maxIts = 100;

%Material/problem parameters
parameters.Q = 1; %or 2 m^3/ms ?
parameters.drho = drho;
%parameters.g = 9.81;
parameters.Kc = 1;%*Ks;

%% define initial conditions (h0, zf0)
dz = parameters.dz;
z0 = -100;
z = (0:dz:11.8)';
n = length(z);


M = csvread("initial_width.csv");

%interpolated values of h at z query points
hq = interp1(M(:,1),M(:,2),z);
hq2 = hq;
tf = (z < 8);
hq(tf) = 1.0177;


%smoother h:
p = polyfit(z, hq, 9);
v = polyval(p, z);
hq(~tf) = v(~tf);
hq = hq - 0.0177;
hq(n) = 0;

z2 = z + z0;
zf = z2(n) + dz/100;

x = zf - z2;

close all
figure; hold on
plot(z2(1:n-1),hq(1:n-1),'o'); 
%plot(z(1:n-1),v(1:n-1),'o-'); 
%M1 = matrixM(x,dz);
[M1,M2] = matrixMreduced(x,dz);
%% 

%Pressure profile that has an analytic soltion
lambda = 2.36;
Pe0 = ones(n,1);
tf = (x > lambda);
Pe0(tf) = 0;

%corresponding h profile according to h = M1Pe0
h0 = M1*Pe0;
%analytic profile
han = analyticSolution(x,lambda);

diff = h0./han;

close all
figure(1); hold on
plot(x,h0,'-o')
plot(x,han,'-k')
%plot(x,diff,'r')
plot(z,Pe0,'r')

Minv = real(han)/real(Pe0);

Pe1 = Minv*hq;

plot(x,Pe1,'b')

%%
tol = 1e-02;
maxit = 499;
% %Pe_gm = gmres(M1s,h2b(1:nz-1),[],tol,maxit); 
hn = (2^0.5);
hq2 = hq;
hq2(n) = hn;
close all

%M = zeros(N,N);
d1 = M1(  1:1+n:n*n); %diagonal k = k
d2 = M1(n+1:1+n:n*n); % k - 1
d3 = M1(  2:1+n:n*n-n); %k + 1

%Jacobi block
J = zeros(n,n);
J(1:1+n:n*n) = d1;
J(n+1:1+n:n*n) = d2;
J(2:1+n:n*n-n) = d3;

%ilu
%J = ichol(M1);
%

figure(1); pc = pcolor(M1)
colorbar
pc.LineStyle='none'

%
res = [];
Pe = gmres(M1,hq2,res,tol,maxit);
PeJ = gmres(M1,hq2,res,tol,maxit,J);
Minv = inv(M1);
Minv2 = Pe/hq2;

Pe3 = Minv2*hq2;
% close all
figure(2); hold on

h2 = M1*Pe;
h3 = M1*Pe1;
plot(z,h2,'s');
plot(z,hq2,'-k');

x1 = 1;
x2 = n;

figure(3); hold on
plot(z(x1:x2),Pe(x1:x2),'-ro');
% plot(z,PeJ,'-g*');

figure(4); pc = pcolor(M1(x1:x2,x1:x2))
colorbar
pc.LineStyle='none'
%%
close all
[h,hmax,zf,t] = dykeModel(parameters,hq,z2); 

%%
% figure(1);
% plot(z,h);
% 
 figure(3);
 plot(1:parameters.nTimeSteps,hmax)
% 
figure(4);
plot(1:parameters.nTimeSteps,fnorm)
%xlabel("Time step")
%ylabel("norm(f)")



