function [LS,V] = constructLSV(Minv,V0,Pe,Kc)
%function to obstain matrix LS and vector V in Pe = LS*h + V

n = length(Pe);

V = Kc*(2^0.5)*Minv(:,n);

Pe0 = Pe;
Pe0(1:n-2) = 1;

% V0(1:n-1) = (M2(1:n-1,1:n)*Pe0);
% V0(n) = M2(n,:)*Pe;
% V0 = Minv*(V0');

LS = Minv;

LS(:,n) = -V0;

