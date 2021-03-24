function dykeModel(parameters, hin, zin, Kc, filename)

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
zf = z(n) + dz/10;

x = zf - z;
x0 = x(1)/xa + 50;

%BC for h at n-1
h(n-1) = Kc*(2^0.5)*(x(n-1))^0.5; % h at n-1 
h(n) = 0;

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
     h = real(h);
     %update zf and Pe
     zf = zf0 + h(n);     
     Pe = elasticPressure(h,LS,V);
     Pe = real(Pe);
     %check if zf > z(n), and add extra node if so
     if zf > z(n) + dz      
         z = vertcat(z, z(n) + dz);%  
         [h,n,Pe] = addNodes(z,h,Pe0,dz,zf,Kc);  
     end  
     if zf < z(n)
        z2 = z(1:n-1);
        h2 = h(1:n-1);
        Pe = Pe(1:n-1);
        n = n-1;
        h2(n-1) = Kc*(2^0.5)*(zf-z2(n-1))^0.5;
        h2(n) = h(n+1);
        h = h2;
        z = z2;
     end
             
     if mod(i,nPlot) == 0
        fprintf('Completed step %d \n', i)
        fprintf('max h is %.3f \n', max(h))
%          fig1 = figure(1); hold on
%          plot(vertcat(z(1:n-1),zf),vertcat(h(1:n-1),0),'o-');  
%          plot(z,Pe,'--')        
%          fig2 = figure(2); hold on
%          plot(t,zf,'-o')
%          xlabel("Time")
%          ylabel("zf")
%          fig3 = figure(3); hold on
%          E = (sum(h.^2))^0.5;
%          plot(t,max(h),'ko-');
         
         file = ['output/' filename '_' num2str(i) '.txt'];
         z_out = vertcat(z(1:n-1),zf);
         h_out = vertcat(h(1:n-1),0);
         M_out = horzcat(z_out,h_out,Pe);
         save(file,'M_out','-ascii');
         
      end

 end

end 




