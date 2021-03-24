function [h,n,Pe] = addNodes(z0,h0,Pe0,dz,zf,Kc)

%function to add extra node to z when zf - z(n) > dz

%increase position vector z by dz
n = length(h0);
% znew = z0(n) + dz;
%z = vertcat(z0(1) - dz,z0);
% z = vertcat(z0,z0(n) + dz);

h = vertcat(h0,h0(n));
h(n) = Kc*(2^0.5)*(zf-z0(n))^0.5;
%h(n-1) = Kc*(2^0.5)*(zf-z0(n-1))^0.5;

Pe = vertcat(Pe0,Pe0(n)); %add extra value of Pe = 0 to the bottom of the dyke
Pe(n) = 0.5*(Pe(n-1) + Pe(n+1));

%add extra point to h, keep the last value of h the same
%h = vertcat(h0, h0(n));
% update n and boundary condition 
 n = n + 1;
 %h(n) = z(n) - z0(n-1)
 %h(n-1) = Kc*(2^0.5)*(zf-z(n-1))^0.5;