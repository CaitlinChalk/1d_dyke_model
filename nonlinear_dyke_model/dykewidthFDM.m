function f = dykewidthFDM(parameters,h,h0,Pe,Pe0)

%function defining the discretised, nonlinear equation for dyke width:
%f(h) = 0
%input - h - size n array of dyke width. h(n) = dzf
%      - h0 - value of h at previous time step
%      - parameters - model input parameters
%      - Pe - elastic pressure (n x 1)
%      - Pe0 - elastic pressure at previous time step
%output - f

n = length(h);
%f = zeros(n-3,1);

dz = parameters.dz;
dt = parameters.dt;
Q = 1; %parameters.Q; %flux inlet BC

n = length(h);

phi2 = (h(2)^3)*( (Pe(2) - Pe(1))/dz - 1);
phi02 = (h0(2)^3)*( (Pe0(2) - Pe0(1))/dz - 1);

%f(1) = h(2) - h0(2) + (dt/(2*dz))*(phi2 + phi02 + 2*Q);

f = h(2:n-2) - h0(2:n-2) - (dt/(16*(dz^2)))*(fdiff(h,Pe,dz,n) + fdiff(h0,Pe0,dz,n)); 

%f(1) = h(2) - h0(2) + (dt/(2*dz))*(phi2 + phi02 + 2*Q);

end


function v = fdiff(h,Pe,dz,n)
    
v = (h(2:n-2) + h(3:n-1)).^3.*(Pe(3:n-1) - Pe(2:n-2) - dz)  - (h(2:n-2) + h(1:n-3)).^3.*(Pe(2:n-2) - Pe(1:n-3) - dz);

end