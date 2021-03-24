function [h,hmax,z,zf] = dykeModel(parameters, hin, zin, zf, Kc)

%solves the governing dyke model equations using Newton's method
%input - parameters - model parameters
%      - hin - initial values of h
%      - z - array of spatial coordinates
%      - nTimeSteps - no. of time steps

%Model parameters
dz = parameters.dz;
xa = parameters.xa;

%initial data for h, t and z
t = 0;
h = hin;
% z = zin;

n = length(h);
z = zin(1:n);
%zf = z(n) + dz/2; %initialise zf
%zf00 = zf; %store initial value of zf
%z0 = min(z);
x = zf - z;
x0 = x(1)/xa + 50;

%BC for h at n-1
h(n-1) = Kc*(2^0.5)*(x(n-1))^0.5; % h at n-1 

%initial pressure profile
Pe0 = zeros(n,1);
[M1,V1] = constructM1V1(x,Pe0,dz,Kc,x0);
h2 = h;
h2(n) = Kc*(2^0.5);
Pe = M1\h2;

%Convergence parameters
nTimeSteps = parameters.nTimeSteps;
nPlot = parameters.nPlot;
dt = parameters.dt;

%output arrays
hmax = zeros(1,nTimeSteps);
figure(1); hold on
%plot(x(1:n-1),h(1:n-1),'k','LineWidth',2)
xlabel("zf - z")
ylabel("h")
figure(2); hold on
plot(0,zf,'k*')

%MAIN LOOP
 for i = 1:nTimeSteps
     t = t + dt;
     %initialise variables
     h0 = h;
     zf0 = zf;
     Pe0 = Pe;
     x = zf0 - z; 
     %Construct M1 and V1 (h = M1Pe + dz*V1)
     [M1,V1] = constructM1V1(x,Pe0,dz,Kc,x0);
     %Calculate the inverse of M1 (?)
     A = decomposition(M1);
     I = eye(n);
     Minv = A\I;
     %Construct LS and V (Pe = LSh + V)
     [LS,V] = constructLSV(Minv,V1,Pe0,Kc);
     %compute A and b and solve Ah = b
     h = solveSystem(parameters, h0, Pe0, LS, V, x(n),x(n-1),Kc);
     %update zf and Pe
     zf = zf0 + h(n);     
     Pe = elasticPressure(h,LS,V);
     %check if zf > z(n), and add extra node if so
     if zf > z(n) + dz      
         z = vertcat(z, z(n) + dz);%  
         [h,n,Pe] = addNodes(z,h,Pe0,dz,zf,Kc);  
     end   
             
     if mod(i,nPlot) == 0
        fprintf('Completed step %d \n', i)
        fprintf('max h is %.3f \n', max(h))
         fig1 = figure(1); hold on
         plot(vertcat(z(1:n-1),zf),vertcat(h(1:n-1),0),'o-');  
         plot(z,Pe,'--')        
         fig2 = figure(2); hold on
         plot(t,zf,'-o')
         xlabel("Time")
         ylabel("zf")
         fig3 = figure(3); hold on
         E = (sum(h.^2))^0.5;
         plot(t,max(h),'ko-');
%          xlabel("Time")
%          ylabel("width")
%          saveas(fig1,['profile' num2str(Kc) '_' num2str(i) '.jpg'])
%          saveas(fig2,['tip' num2str(Kc) '_' num2str(i) '.jpg'])
%          saveas(fig3,['hmax' num2str(Kc) '_' num2str(i) '.jpg'])
      end
%      if mod(i,100000) == 0
%         figure(4); hold on
%         plot(z,vertcat(h(1:n-1),0),'k-');
%      end
 end

end 




