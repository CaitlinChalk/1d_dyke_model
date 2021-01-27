function dzf = tipIncrement2(parameters,h,h0)

%function to calculate the tip increment using volume conservation across
%the whole dyke
%input - parameters - problem parameters
%      - h - dyke width
%      - h0 - dyke width at previous time step 

%output - dzf - dyke front increment between current and previous time step

Q = parameters.Q;
dt = parameters.dt;
n = length(h);

h2 = h0(1:n-1) - h(1:n-1);

h2int = -trapz(h2);

dzf = (2/h(n-1))*(Q*dt + h2int);


