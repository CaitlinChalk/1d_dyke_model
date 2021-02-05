function [h,hmax,zf,t] = dykeModel(parameters, hin, zin)

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
zf = z(n) + dz/10; %initialise zf
z0 = min(z);
x = zf - z;
zf2 = zf; %save separate zf to compare two different calculation methods for h(n)/dzf
zf3 = zf;

%BC for h at n-1
h(n-1) = Kc*(2^0.5)*(x(n-1))^0.5; % h at n-1 

%initial pressure and Minv
warning('off','MATLAB:gmres:tooBigTolerance')
[M1,M2] = matrixMreduced(x,dz);
tol = 1e-6;
maxit = n-1;
h2 = h;
h2(n) = Kc*(2^0.5);
Pe = gmresnomsg(M1,h2,[],tol,maxit);
%Pe = elasticPressure(parameters, h, Pe0, M1,M2); %initial guess for pressure
Minv = Pe/h2;


%Convergence parameters
tol = parameters.tol;
maxIts = parameters.maxIts;
nTimeSteps = parameters.nTimeSteps;
dt = parameters.dt;

%output arrays
hmax = zeros(1,nTimeSteps);
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
     [h,Pe,f] = newtonAlgorithm(parameters, h0, x(n), Pe0, Minv, M2, @dykewidthFDM,@fdJacobian,tol, maxIts);     
     Minv = updateMinv(M1,h,Kc);
     %h(1) = (1/( (Pe(2)-Pe(1))/dz + 1))^(1/3);
     h(n) = tipIncrement(parameters,h,h0,Pe,Pe0,x(n));    
     %h(n) = tipIncrement2(parameters,h,h0);
     zf = zf0 + h(n);
   %  zf2 = zf2 + hn2;
%      zf3 = zf3 - hn2;
     %update BC for h at n-1    
     hmax(i) = max(h(1:n-1));
     if mod(i,1000) == 0
        fprintf('Completed step %d \n', i)
        fprintf('norm(f) is %.4f \n', norm(f))
        fprintf('max h is %.3f \n', hmax(i))
        figure(1); hold on
        plot(x,vertcat(h(1:n-1),0),'o-');
        plot(x(1:n-1),Pe(1:n-1),'-')        
        figure(2); hold on
        plot(t,zf,'-o')
      %  plot(dt*i,zf2,'-.^')
%         plot(dt*i,zf3,'-.*')
        xlabel("Time")
        ylabel("zf")
        figure(3); hold on
        plot(t,hmax(i),'o-');
        xlabel("Time")
        ylabel("hmax")
     end
     if mod(i,100000) == 0
        zn = convertToZ(z0,zf,n);
        figure(4); hold on
        plot(zn,vertcat(h(1:n-1),0),'k-');
     end
 end

end 

function x = gmresnomsg(varargin)
    [x,~] = gmres(varargin{:});
end




