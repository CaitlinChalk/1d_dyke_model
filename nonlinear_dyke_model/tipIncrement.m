function dzf = tipIncrement(parameters,h,h0,Pe,Pe0,z,zf0)

%function to calculate the change in tip position betwen zf|t and zf|t+1
%input - parameters - problem parameters
%      - h - dyke width
%      - h0 - dyke width at previous time step
%      - z - height
%      - zf0 - dyke tip position at previous time step

%output - dzf - change in tip position 

n = length(h);

dt = parameters.dt;
hs = parameters.hs;
dz = parameters.dz;

hn2 = h0(n-2)/h(n-2);

dzf = z(n) - zf0 + (dt/(2*hs^2*dz))*(hn2*(Pe0(n-1) - Pe0(n-2) + zf0 - z(n) - dz) + Pe(n-1) - Pe(n-2) - dz );


