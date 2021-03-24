function [M1,V1] = constructM1V1(x,Pe,dz,Kc,x0)

%function to calculate the matrix M1 and vector V1

%% M1

% x coordinate system

n = length(x);
M1 = zeros(n);
V1 = zeros(n,1);

for i = 1:n-1
    for j = 1:n
      
%         A(i,j) = getA(x(i),x(j),dz);
%         B(i,j) = getB(x(i),x(j),dz);
%         f(i,j) = getF(x(i),x(j));
%         g(i,j) = getG(x(i),x(j));
        
        if j == 1                       
            M1(i,j) = getA(x(i),x(j),x(j+1),dz) + x0/(x0-x(1))*(getF(x(i),x(1)) - getF(x(i),x0)) -  (getG(x(i),x(1)) - getG(x(i),x0))/(x0-x(1)); 
        elseif j == n-1
            M1(i,j) = getB(x(i),x(j-1),x(j),dz) - getG(x(i),x(j))/x(j);
        elseif j == n
            M1(i,j) = -getF(x(i),x(j-1)) + getG(x(i),x(j-1))/x(j-1); 
        else
            M1(i,j) = getA(x(i),x(j),x(j+1),dz) + getB(x(i),x(j-1),x(j),dz);
        end        
    end
    V1(i)=Kc*((2.0/x(i))^0.5)*pi/2.0+(Pe(n)-Pe(n-1))*(-getG(x(i),x(n-1))/x(n-1)^2.0+getF(x(i),x(n-1))/x(n-1));
end

M1(n,1) = (1 - x(1)/dz)*(getFtip(x(2)) - getFtip(x(1))) + (getGtip(x(2)) - getGtip(x(1)))/dz;
M1(n,1) = M1(n,1) + x0/(x0-x(1))*(getF(x(n),x(1)) - getF(x(n),x0)) - (getG(x(n),x(1)) - getG(x(n),x0))/(x0-x(1));
% M1(n,n-1) = -getGtip(x(n-1))/x(n-1) + (x(n-2)/dz)*(getFtip(x(n-1)) - getFtip(x(n-2))) - (getGtip(x(n-1)) - getGtip(x(n-2)))/dz;
% M1(n,n) = -getFtip(x(n-1)) + getGtip(x(n-1))/x(n-1);
M1(n,n-1) = (x(n-2)./dz).*(getFtip(x(n-1)) - getFtip(x(n-2))) - (getGtip(x(n-1)) - getGtip(x(n-2)))./dz - getGtip(x(n-1))/x(n-1);
M1(n,n) = -getFtip(x(n-1)) + getGtip(x(n-1))/x(n-1);

M1(n,2:n-2) = (1 - x(2:n-2)./dz).*(getFtip(x(3:n-1)) - getFtip(x(2:n-2))) ...
                + (getGtip(x(3:n-1)) - getGtip(x(2:n-2)))./dz ...
                + (x(1:n-3)./dz).*(getFtip(x(2:n-2)) - getFtip(x(1:n-3))) - (getGtip(x(2:n-2)) - getGtip(x(1:n-3)))./dz;

            
%V1(i)=V1(i)+pe(n-1)/(3.0_dp*dsqrt(zf-z(n-1)))+2.0_dp*pe(n)/(3.0_dp*dsqrt(zf-z(n-1)))
V1(n)=V1(n)+5.0/(6.0*(x(n-1)^0.5))*(Pe(n-1)-Pe(n));
            
            
M1 = (1/pi).*M1;
V1 = (1/pi).*V1;

M1(n,:) = 2.*M1(n,:);
V1(n) = 2.*V1(n);
%z coordinate system




end

%functions required to construct matrices M1 and M2

function Aij = getA(i,j,jp,dz)
    Aij = (1 - (j/dz))*( getF(i,jp) -  getF(i,j) ) + (getG(i,jp) - getG(i,j))/dz;
end

function Bij = getB(i,j,jp,dz)
    Bij = (j/dz)*( getF(i,jp) -  getF(i,j) ) - (getG(i,jp) - getG(i,j))/dz ; 
end

function fij = getF(i,j)
    fij = (i-j).*log( abs( (i.^0.5 + j.^0.5)./(i.^0.5 - j.^0.5) ) ) - 2.*(i.*j).^0.5; 
    if i == j
        fij =  - 2.*(i.*j).^0.5;
    end
end

function gij = getG(i,j)
    gij = 0.5.*(i.^2-j.^2).*log( abs( (i.^0.5 + j.^0.5)./(i.^0.5 - j.^0.5) ) ) - (i.^0.5).*(j.^(3/2))./3 - (i.^(3/2)).*(j.^0.5); 
    if i == j
        gij =  -(i.^0.5).*(j.^(3/2))./3 - (i.^(3/2)).*(j.^0.5);
    end
end

function ftip = getFtip(j)
    ftip = -2*(j.^0.5);
end

function gtip = getGtip(j)
    gtip = -(2/3)*j.^(3/2);
end


function Aj = getAn(j,jp,dz)
    Aj = (1 - j./dz).*((-1./(jp.^0.5))+(1./(j.^0.5)) ) +(jp.^0.5 - j.^0.5)./dz;
end

function Bj = getBn(j,jp,dz)
    Bj = (j./dz).*((-1./(jp.^0.5))+(1./(j.^0.5)) ) -(jp.^0.5 - j.^0.5)./dz;
end




