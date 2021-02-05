function [M1,M2] = matrixMreduced(x,dz)

%function to calculate the matrices M1 and M2

%% M1

% x coordinate system

n = length(x);
M1 = zeros(n);

for i = 1:n-1
    for j = 1:n


        if i ~= j 
            k = klogterm(x(i),x(j));
        else
            k = 0;
        end
        
        if j == 1           
            M1(i,j) = k; 
        else
            M1(i,j) = 2*k;
        end

    end
end


M1(n,2:n-1) = 4./(x(2:n-1).^0.5)';
M1(n,1) = 2/(x(1)^0.5);
M1(n,n) = 0;

M1 = (dz/(2*pi))*M1;


%z coordinate system

% n = length(z);
% M1 = zeros(n);
% dz = z(n) - z(n-1);
% 
% for i = 1:n-1
%     for j = 1:n
%         if i ~= j            
%             M1(i,j) = 2*klogterm(z(i),z(j),zf);    
%         end
%     end
% end
% 
% M1(n,:) = (4./((zf - z).^0.5))';
% M1(:,[1,n]) = 0.5*M1(:,[1,n]);
% M1 = (dz/(2*pi))*M1;



%% M2

M2 = zeros(n-1);

for i = 1:n-1
    for j = 1:n

        km = kdiff(x(i),x(j) + dz/2);
        kp = kdiff(x(i),x(j) - dz/2);

        if i ~= j 
            k = kdiff(x(i),x(j));
        else
            k = 0;
        end
        
        if j == 1           
            M2(i,j) = k - kp; 
        elseif j == n
            M2(i,j) = km;
        else
            M2(i,j) = km + 2*k - kp;
        end

    end
end

M2(n,2:n-1) = (-2./((x(2:n-1)+dz/2).^1.5) - 4./(x(2:n-1).^1.5) + 2./((x(2:n-1)-dz/2).^1.5) )';
M2(n,1) = -4/(x(1)^1.5) + 2/((x(1)-dz/2.)^1.5);
M2(n,n) = -2/((x(n)+dz/2.)^1.5);


M2 = (dz/(4*pi))*M2;


% M2 = zeros(n-1);
% 
% for i = 1:n-1
%     for j = 1:n
%         M2(i,j) = 2/((zf - z(i))*(zf - z(j)));
%     end
% end
% 
% M2(n,:) = -2./((zf - z).^(3/2));
% M2(:,[1,n]) = 0.5*M2(:,[1,n]);
% M2 = (dz/(2*pi))*M2;

end
function kij = klogterm(a,b)
    kij = log( abs( (a^0.5 + b^0.5)/(a^0.5 - b^0.5) ) ); 
end

function kij2 = kdiff(a,b)
    kij2 = 1/((a^0.5)*(b^0.5));
end
% function kij = klogterm(i,j,zf)
%     kij = log( abs( ((zf - i)^0.5 + (zf - j)^0.5)/((zf - i)^0.5 - (zf - j)^0.5) ) ); 
% end