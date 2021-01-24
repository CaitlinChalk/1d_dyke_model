function [z,h,Pe] = addNodes(z0,h0,Pe0,dz)

%function to add extra node to z when zf - z(n) > dz
n = length(z0);
znew = z0(n) + dz;
z = vertcat(z0,znew);

%interpolate h0 and Pe0
%h2 = interp1(z0(1:n-1),h0(1:n-1),z(1:n));
%Pe = interp1(z0,Pe0,z);

h = vertcat(h0, h0(n));
Pe = vertcat(Pe0, Pe0(n));
%Pe = ones(n+1,1);

%boundary condition
h(1) = h0(1);

%to do - case when zf - z(n) > 2dz etc