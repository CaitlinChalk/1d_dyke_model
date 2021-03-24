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

dz = 0.0472/2;

%define input parameters to pass to model
parameters = struct();


%Material/problem parameters
parameters.Q = Q; %or 2 m^3/ms ?
parameters.drho = drho;
%parameters.g = 9.81;

%dimensionless numbers
dza = dz;
ha = (3*mu*Q/(2*drho*g))^(1/3);
ca = Q/(2*ha);
%xa = ((3*mu*Q/2)^(1/6))*((G/(1-nu))^0.5)*(drho*g)^(2/3);
xa=(G*ha/((1.0-nu)*(drho)*g))^0.5;
ta = xa/ca;
dta = (20/500)*(dza);%dza = dz/za

%numerical parameters
parameters.dt = dta;
parameters.dz = dz;
parameters.nTimeSteps = 2000000;
parameters.nPlot = 10;%round(0.1/parameters.dt);


parameters.ha = ha;
parameters.ca = ca;
parameters.xa = xa;
parameters.ta = ta;




% define initial conditions (h0, zf0)

z0 = -100;
z = (0:dza:11.8)';
n = length(z);

%initial condition 1 - classic dyke shape
%read in h profile from file (extrapolated from taisne & juapart 2009,fig
%C1
M = csvread("initial_width.csv");

%interpolated values of h at z query points
hq = interp1(M(:,1),M(:,2),z);
%manually smooth profile
tf = (z < 8);
hq(tf) = 1.0177;

%smoother h:
p = polyfit(z, hq, 9);
v = polyval(p, z);
hq(~tf) = v(~tf);
h1 = hq - 0.0177;
h1(n) = 0;

dz0 = dza/100;
z2 = z + z0;
zf = z2(n) + dz0;

z_full = (z0:dz:z0+59)';
%initial condition 2 - blunt dyke

h2_kc1 = ones(n,1); %for Kc = 1
h2_kc2 = ones(n,1); %for Kc = 2

check1 = 1.*(2^0.5).*((zf - z2).^0.5); %Kc = 1;
check2 = 2.*(2^0.5).*((zf - z2).^0.5); %Kc = 2;

h2_kc1(check1 < 1) = check1(check1 < 1);
h2_kc1(n) = 0;

h2_kc2(check2 < 1) = check2(check2 < 1);
h2_kc2(n) = 0;
%coordinate reference system that moves with the front tip position
x = zf - z2;

close all
figure(1); hold on
plot(z2 - z0,h1,'o-'); 
plot(z2 - z0,h2_kc1,'o-'); 
xlabel("$z - z_0$","Interpreter","Latex")
ylabel("$h$","Interpreter","Latex")
title("Initital dyke width profile","Interpreter","Latex")
legend(["Initial condition 1","Initial condition 2"],"Interpreter","Latex","Location","NorthWest")


%% Initial pressure profile

%x = b/A => xA = b
%x = A\b => Ax = b
%h = M1*Pe => M1 = h/Pe, Pe = M1\h;
%Minv*h = Pe => Minv = Pe/h

Pe0 = zeros(n,1);
[M1,V1] = constructM1V1(x,Pe0,dz,1,x(1)/xa + 50);

hn = (2^0.5);
hq2 = hq; %hq
hq2(n) = hn;

%numerical solution using \ operator
Pe_numeric = M1\hq2;

%numerical solution using inv
Minv = inv(M1);
Pe_numeric2 = Minv*hq2;

%alternative way to compute inv(M1)
A = decomposition(M1);
I = eye(n);
Minv2 = A\I;

%corresponding h profile (back calculation for validation purposes)
h_numeric = M1*Pe_numeric; 

close all
figure(2); hold on
plot(x,h_numeric,'k','LineWidth',2)
plot(x,Pe_numeric,'r-o')

xlabel("$z_f - z$","Interpreter","Latex")
ylabel("$P_e$","Interpreter","Latex")
title("Initital pressure profile","Interpreter","Latex")



%% run model
close all

%initial condition 1, Kc = 1
[h_final,hmax,z,zf] = dykeModel(parameters,h1,z_full,zf,1); 
% ic1_kc1 = {parameters,h1,z2,1};
% 
% %initial condition 2, Kc = 1
% % [h_final,hmax,z,zf] = dykeModel(parameters,h2_kc1,z2,1); 
% ic2_kc1 = {parameters,h2_kc1,z2,1};
% 
% % %initial condition 2, Kc = 2
% % [h_final,hmax,z,zf] = dykeModel(parameters,h2_kc2,z2,2); 
% ic2_kc2 = {parameters,h2_kc2,z2,2};
%% parallel processes
% func = @dykeModel;
% arguments = {parameters,h1,z2,1 ; parameters,h2_kc1,z2,1 ; parameters,h2_kc2,z2,2};
% h_final = cell(1,3);
% h_max = cell(1,3);
% z_final = cell(1,3);
% zf = cell(1,3);
% 
% %
% parfor i = 1:2
%     [hi,hmi,zi,zfi] = func(arguments{i,:});
% 
%   %  [h_final(i,1:4)] = func(arguments{i,:})
%      h_final{i} = hi;
%      h_max{i} = hmi;
%      z_final{i} = zi;
%      zf{i} = zfi;
% end



