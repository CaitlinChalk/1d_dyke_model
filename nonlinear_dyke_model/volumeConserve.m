function dzf = volumeConserve(h,h0,xn,dz,dt,Q,Kc)


n1 = length(h);
n0 = n1 - 1;

%dzf = (1/(0.5*h0(n0-2) + Kc*(dz^0.5)))*(Q*dt + 0.5*h0(1)*dz + dz*sum(h0(2:n0-3)) + 0.5*h0(n0-2)*(dz + xn) ...
%       -(0.5*h(1)*dz + dz*sum(h(2:n1-3)) +  0.5*h(n1-2)*(dz + xn) ));
   
   
dzf = (2/(dz + xn))*(Q*dt + 0.5*h0(1)*dz + dz*sum(h0(2:n0-3)) + 0.5*h0(n0-2)*(dz + xn)...
     - ( 0.5*h(1)*dz + dz*sum(h(2:n1-3)) + h(n1)*(0.5*h0(n0-2) + Kc*(dz^0.5))));