function [h,hmax,fnorm,t] = dykeModel(parameters, hin, zin)

%solves the governing dyke model equations using Newton's method
%input - parameters - model parameters
%      - hin - initial values of h
%      - z - array of spatial coordinates
%      - nTimeSteps - no. of time steps

%Model parameters
Kc = parameters.Kc;
dz = parameters.dz;

%initial data for h, t and z
t = 0;
h = hin;
z = zin;

n = length(h);
zf = z(n) + dz/100; %initialise zf
x = zf - z;

%update BC for h at n-1
h(n-1) = Kc*(2^0.5)*(x(n-1))^0.5; % h at n-1 

%initial pressure and Minv
Pe0 = zeros(n,1);
[M1,M2] = matrixM(x,dz);
Pe = elasticPressure(parameters, h, Pe0, M1,M2); %initial guess for pressure

%Convergence parameters
tol = parameters.tol;
maxIts = parameters.maxIts;
nTimeSteps = parameters.nTimeSteps;
dt = parameters.dt;

%output arrays
hmax = zeros(1,nTimeSteps);
fnorm = zeros(1,nTimeSteps);
figure(1); hold on
plot(x(1:n-1),h(1:n-1),'k','LineWidth',2)
xlabel("zf - z")
ylabel("h/Pe")
figure(2); hold on
plot(0,zf,'k*')
 for i = 1:nTimeSteps
     t = t + dt;
     h0 = h;
     zf0 = zf;
     Pe0 = Pe;
     [h,Pe,f] = newtonAlgorithm(parameters, h0, Pe0, M1, M2, @dykewidthFDM,@fdJacobian,tol, maxIts);      
     h(n) = tipIncrement2(parameters,x(n-1),zf0,h,h0);    
     zf = zf0 + h(n);
     %update BC for h at n-1
     h(n-1) = Kc*(2^0.5)*(x(n-1))^0.5; % h at n-1   
     hmax(i) = max(h(1:n-1));
     fnorm(i) = norm(f);
     if mod(i,500) == 0
        fprintf('Completed step %d \n', i)
        fprintf('norm(f) is %.4f \n', norm(f))
        fprintf('max h is %.3f \n', hmax(i))
        figure(1)
        plot(x,vertcat(h(1:n-1),0),'o-');
        plot(x(1:n-1),Pe(1:n-1),'-')        
        figure(2); hold on
        plot(dt*i,zf,'-o')
%         plot(i,zf2,'-.^')
        xlabel("Time")
         ylabel("zf")
         figure(3); hold on
        plot(dt*i,hmax(i),'o-');
        xlabel("Time")
         ylabel("hmax")
     end
     
 end




