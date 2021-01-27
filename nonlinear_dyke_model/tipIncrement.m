function dzf = tipIncrement(parameters,h,h0,Pe,Pe0,xn)

%function to calculate the change in tip position using volume conservation
%in the near-tip domain only

%input - parameters - problem parameters
%      - h - dyke width
%      - h0 - dyke width at previous time step
%      - Pe - dyke pressure
%      - Pe0 - dyke pressure at previous time step
%      - xn - nth value of position vector x

%output - dzf - change in tip position 

n = length(h);

dt = parameters.dt;
dz = parameters.dz;

phi = h(n-2)^3*(-(Pe(n-1) - Pe(n-3))/(2*dz) + 1);
phi0 = h0(n-2)^3*(-(Pe0(n-1) - Pe0(n-3))/(2*dz) + 1);

dzf =  (xn/h0(n-2))*(h0(n-2) - h(n-2)) + (dt/(2*h0(n-2)))*(phi0 + phi); 


