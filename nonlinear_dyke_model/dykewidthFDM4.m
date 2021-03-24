function [A,b] = dykewidthFDM4(parameters,h0,Pe0,LS,V,xn,xn1,Kc)

%function defining the discretised, nonlinear equation for dyke width:
%f(h) = 0
%input - h - size n array of dyke width. h(n) = dzf
%      - h0 - value of h at previous time step
%      - parameters - model input parameters
%      - Pe - elastic pressure (n x 1)
%      - Pe0 - elastic pressure at previous time step
%output - f

n = length(h0);
A = zeros(n,n);
b = zeros(n,1);

dz = parameters.dz;
dt = parameters.dt;
Q = parameters.Q; %flux inlet BC

%from thesis eq 5, constant flux, 

% A(1,:) = ((h0(1)^3)*((LS(2,:)-LS(1,:))/dz));
% A(1,1) = A(1,1) + (3*h0(1)^2)*((Pe0(2)-Pe0(1))/dz - 1);
%   
% b(1) = -(Q - (h0(1)^3)*(3*(Pe0(2)-Pe0(1))/dz -2 - (V(2)-V(1))/dz));

%from fortran code, constant flux 

A(1,1)=3.0*h0(1)^2.0*((Pe0(2)-Pe0(1))/dz-1.0);

A(1,1:n)=A(1,1:n)+(h0(1)^3.0).*((LS(2,1:n)-LS(1,1:n))./dz);

b(1)=-Q+h0(1)^3.0*(3.0*(Pe0(2)-Pe0(1))/dz-2.0-(V(2)-V(1))/dz);


%from thesis
% for i = 2:n-2
%     A(i,:) = (-fdiff(LS(i+1,:),LS(i-1,:),dz)*((3/2)*dt*(h0(i)^2)*fdiff(h0(i+1),h0(i-1),dz)) ...
%         - fdiff2(LS(i,:),LS(i+1,:),LS(i-1,:),dz).*(dt/2).*h0(i)^3);
%     A(i,i) = A(i,i) ...
%         + (1 - ((3*dt)/2)*(h0(i)^2)*fdiff2(Pe0(i),Pe0(i+1),Pe0(i-1),dz) - 3*dt*h0(i)*fdiff(h0(i+1),h0(i-1),dz)*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1));
%     A(i-1,i) = A(i-1,i) + ((3*dt/4*dz)*(h0(i)^2)*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1));
%     A(i+1,i) = (-(3*dt/4*dz)*h0(i)^2*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1));
%             
%     b(i) = h0(i) + dt*(h0(i)^2)*fdiff(h0(i+1),h0(i-1),dz)*(-3*fdiff(Pe0(i+1),Pe0(i-1),dz) + 3/2) ...
%     - dt*(h0(i)^3)*fdiff2(Pe0(i),Pe0(i+1),Pe0(i-1),dz) ...
%     + (3*dt/2)*(h0(i)^2)*fdiff(h0(i+1),h0(i-1),dz)*fdiff(V(i+1),V(i-1),dz) ...
%     + (dt/2)*(h0(i)^3)*fdiff2(V(i),V(i+1),V(i-1),dz);  
%  
% end

%from fortran code

for i=2:n-2
      A(i,:)=A(i,:)-(LS(i+1,:)-LS(i-1,:))*(3.0*dt*h0(i)^2.0*fdiff(h0(i+1),h0(i-1),dz)/(4.0*dz)) ...
                     -(LS(i+1,:)-2.0*LS(i,:)+LS(i-1,:))*dt*h0(i)^3.0/(2.0*dz^2.0);

   A(i,i-1)=A(i,i-1)+3.0*dt*h0(i)^2.0*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1.0)/(4.0*dz);
   A(i,i)=A(i,i)+1.0-3.0*dt*h0(i)^2.0*fdiff2(Pe0(i),Pe0(i+1),Pe0(i-1),dz)/2.0-3.0*dt*h0(i)*fdiff(h0(i+1),h0(i-1),dz)*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1.0);
   A(i,i+1)=A(i,i+1)-3.0*dt*h0(i)^2.0*(fdiff(Pe0(i+1),Pe0(i-1),dz)-1.0)/(4.0*dz);

   b(i)=h0(i)+dt*h0(i)^2.0*fdiff(h0(i+1),h0(i-1),dz)*(-3.0*fdiff(Pe0(i+1),Pe0(i-1),dz)+1.5) ...
            -dt*h0(i)^3.0*fdiff2(Pe0(i),Pe0(i+1),Pe0(i-1),dz) ...
            +1.5*dt*h0(i)^2.0*fdiff(h0(i+1),h0(i-1),dz)*(V(i+1)-V(i-1))/(2.0*dz) ...
            +dt*h0(i)^3.0*(V(i+1)-2.0*V(i)+V(i-1))/(2.0*dz^2.0);
end

A(n-1,n-1) = 1;
A(n-1,n) = -(Kc*(2^0.5)/(2*((xn1)^0.5)));
b(n-1) = h0(n-1);

phi = -(h0(n-2)^3)*(fdiff(Pe0(n-1),Pe0(n-3),dz)-1);

%thesis eqn 20 (volume conservation at the tip)

% A(n,:) = ((h0(n-2)^3)*(fdiff(LS(n-1,:),LS(n-3,:),dz))/2);
% A(n,n-2) = A(n,n-2) + xn/(2*dt) + (3*(h0(n-2)^2)/2)*(fdiff(Pe0(n-1),Pe0(n-3),dz) - 1);
% A(n,n) = A(n,n) + (1/dt)*(Kc*(dz^0.5) + (h0(n-2)/2));
% b(n) = h0(n-2)*(xn/(2*dt)) - phi/2 + ((h0(n-2)^3)/2)*fdiff(Pe0(n-1),Pe0(n-3),dz) ...
%        - ((h0(n-2)^3)/2)*fdiff(V(n-1),V(n-3),dz);

%thesis eqn 23 (overall volume conservation)

% A(n,1) = dz/2;
% A(n,2:n-3) = dz;
% A(n,n-2) = (dz + xn)/2;
% A(n,n) = h0(n-2)/2 + Kc*(dz^0.5);
% 
% b(n) = Q*dt + h0(1)*dz/2 + dz*sum(h0(2:n-3)) + h0(n-2)*(dz + xn)/2;
   

%from fortran code (flow condition at the tip)



%from fortran code (global volume conservation)

A(n,1)=A(n,1)+0.5*dz;
b(n)=b(n)+0.5*dz*h0(1);

A(n,2:n-3)=A(n,2:n-3)+dz;
b(n)=b(n)+dz*sum(h0(2:n-3));

A(n,n-2)=A(n,n-2)+0.5*(dz+xn);
A(n,n)=A(n,n)+0.5*h0(n-2)+Kc*(dz^0.5);
b(n)=b(n)+0.5*(dz+xn)*h0(n-2)  + Q*dt;   
   
end

function v = fdiff(xp,xm,dz)
    
v = (xp - xm)./(2*dz);

end


function v2 = fdiff2(x,xp,xm,dz)
    
v2 = (xp - 2.*x + xm)./(dz^2);

end