function dzf = tipIncrement2(parameters,z,h,h0)

%function to calculate the tip increment using volume conservation across
%the whole dyke

Q = parameters.Q;
dt = parameters.dt;
n = length(h);

h2 = h0(1:n-1) - h(1:n-1);

h2int = trapz(h2);

dzf = (2/((1-z(n-1))*h(n-1)))*(Q*dt + h2int);


