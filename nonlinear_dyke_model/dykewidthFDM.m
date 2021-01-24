function f = dykewidthFDM(parameters,h,h0,Pe,Pe0)

%function defining the discretised, nonlinear equation for dyke width:
%f(h) = 0
%input - h - size n array of dyke width. h(n) = dzf
%      - h0 - value of h at previous time step
%      - parameters - model input parameters
%      - Pe - elastic pressure (n x 1)
%      - Pe0 - elastic pressure at previous time step
%output - f

dz = parameters.dz;
dt = parameters.dt;

n = length(h);


f = h(2:n-2) - h0(2:n-2) - (dt/(16*(dz^2)))*(fdiff(h,Pe,dz,n) + fdiff(h0,Pe0,dz,n)); 


end


function v = fdiff(h,Pe,dz,n)
    
v = (h(2:n-2) + h(3:n-1)).^3.*(Pe(3:n-1) - Pe(2:n-2) - dz)  - (h(2:n-2) + h(1:n-3)).^3.*(Pe(2:n-2) - Pe(1:n-3) - dz);

end