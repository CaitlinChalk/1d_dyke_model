function [M1,M2] = matrixM(z,zf)

%function to calculate the matrices M1 and M2

%% M1
n = length(z);
M1 = zeros(n);
dz = z(n) - z(n-1);

for i = 1:n-1
    for j = 1:n
        if i ~= j            
            M1(i,j) = 2*klogterm(z(i),z(j),zf);    
        end
    end
end

M1(n,:) = (4./((zf - z).^0.5))';
M1(:,[1,n]) = 0.5*M1(:,[1,n]);
M1 = (dz/(2*pi))*M1;



%% M2
M2 = zeros(n-1);

for i = 1:n-1
    for j = 1:n
        M2(i,j) = 2/((zf - z(i))*(zf - z(j)));
    end
end

M2(n,:) = -2./((zf - z).^(3/2));
M2(:,[1,n]) = 0.5*M2(:,[1,n]);
M2 = (dz/(2*pi))*M2;

end

function kij = klogterm(i,j,zf)
    kij = log( abs( ((zf - i)^0.5 + (zf - j)^0.5)/((zf - i)^0.5 - (zf - j)^0.5) ) ); 
end