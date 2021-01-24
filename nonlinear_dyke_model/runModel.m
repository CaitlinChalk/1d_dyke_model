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
Q = 2;
Ks = ((3*mu*Q/2)^(1/6))*((G/(1-nu))^(1/2))*((drho*g)^(2/3)); %= 9.7673 10^6
A = 0.1; %scaling factor for Kc: Kc = A*Ks

%define input parameters to pass to model
parameters = struct();

%numerical parameters
parameters.dt = 0.000005;
parameters.dz = 0.0472; %0.0236
parameters.nTimeSteps = 1000000;
parameters.tol = 1.e-6;
parameters.maxIts = 100;

%Material/problem parameters
parameters.Q = 2; %or 2 m^3/ms ?
parameters.drho = drho;
%parameters.g = 9.81;
parameters.hs = ((3*mu*Q)/(2*drho*g))^(1/3);
parameters.alpha = parameters.dt/(4*parameters.hs*(parameters.dz^2));
parameters.Kc = 1;%*Ks;

%% define initial conditions (h0, zf0)
dz = parameters.dz;
z0 = 0;
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
hq(n) = 0;

close all
figure; hold on
plot(z(1:n-1),hq(1:n-1),'o'); 
%plot(z(1:n-1),v(1:n-1),'o-'); 

% [M1,M2] = matrixM(z,z(n) + dz/100);
% 
% tol = 1e-2;
% maxit = 499;
% %Pe_gm = gmres(M1s,h2b(1:nz-1),[],tol,maxit); 
% hn = 2*(2^0.5);
% hq2 = hq;
% hq2(n) = hn;
% Pe = gmres(M1,hq2,[],tol,maxit);
% figure
% plot(z,Pe,'-o')
%%
close all
[h,hmax,fnorm,t] = dykeModel(parameters,hq,z); 

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



