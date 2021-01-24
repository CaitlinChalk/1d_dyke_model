%dyke propagation profiles - Taisne & Jaupart 2009
clear; close all
%Case 1 - unit overpressure

nz = 1000; %no. of grid points

%analytical solution:
dz = 0.00236;
zf = -47 + dz;
lambda = 2.36;
z0 = zf - (dz*nz);

%height
z = linspace(zf - dz,zf - dz - 10*lambda,nz);
%% analytical solution presented in Taisne & Jaupart (doesn't work)
fz = ( (zf - lambda - z).^0.5 + lambda^0.5 )./( (zf - lambda - z).^0.5 - lambda^0.5 );
h = ( 2*(lambda*(zf-z)).^0.5 - (zf - lambda - z).*log(abs(fz)) )/pi;
plot(z,h);

%% wolfram alpha solution to integral (seems to work)
fz = ( ((zf - z).^0.5 + lambda^0.5)./((zf - z).^0.5 - lambda^0.5) );
ft = -2*(zf - z).*atanh( (lambda./(zf-z)).^0.5 );
h0 = ft + lambda*log(abs(fz)) + 2*(lambda*(zf-z)).^0.5;
close all
plot(z,h0);

%% Numerical approximation 1 - trapezium rule with h = dz

%unit pressure vector
Pe = ones(nz,1);

%construct matrix
M1 = zeros(nz);
zr = flip(z);
for i = 1:nz-1
    for j = 1:nz
        if i ~= j            
            M1(i,j) = 2*log(abs( ((zf - zr(i))^0.5 + (zf - zr(j))^0.5)/((zf - zr(i))^0.5 - (zf - zr(j))^0.5)) );    
        end
    end
end

M1(nz,:) = (4./((zf - zr).^0.5))';

%% Numerical approximation 2 - trapezium rule with h = dz/2
%construct matrix
M2 = zeros(nz);
zr = flip(z);

for j = 2:nz-1
    z0 = zr(j);
    z1 = zr(j) - dz/2;
    z2 = zr(j) + dz/2;
    for i = 1:nz-1
        if i ~= j
            k0 = log(abs( ((zf - zr(i))^0.5 + (zf - z0)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z0)^0.5)) );
        else
            k0 = 0;
        end
        k1 = log(abs( ((zf - zr(i))^0.5 + (zf - z1)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z1)^0.5)) );
        k2 = log(abs( ((zf - zr(i))^0.5 + (zf - z2)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z2)^0.5)) );
        M2(i,j) = 2*k0 + k1 + k2; 
    end
end

%first column
z0 = zr(1);
z2 = zr(1) + dz/2;
for i = 2:nz-1
    k0 = log(abs( ((zf - zr(i))^0.5 + (zf - z0)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z0)^0.5)) );    
    k2 = log(abs( ((zf - zr(i))^0.5 + (zf - z2)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z2)^0.5)) );
    M2(i,1) = k0 + k2;
end

%last column
z1 = zr(nz) - dz/2;
for i = 1:nz-1
    k1 = log(abs( ((zf - zr(i))^0.5 + (zf - z1)^0.5)/ ( (zf - zr(i))^0.5 - (zf - z1)^0.5)) );
    M2(i,nz) = k1;
end

%last row

i = nz;
for j = 1:nz
    z0 = zr(j);
    z1 = zr(j) - dz/2;
    z2 = zr(j) + dz/2;
    L0 = 1/((zf - z0)^0.5) ;
    L1 = 1/((zf - z1)^0.5) ;
    L2 = 1/((zf - z2)^0.5) ;  
    if j == 1
        M2(i,j) = 2*(L0 + L1);
    elseif j == nz
        M2(i,j) = L1;
    else
        M2(i,j) = 2*(2*L0 + L1 + L2); 
    end
end

%% plot results

M1b = M1.*(dz/(2*pi));
M2b = M2.*(dz/(4*pi));

Pe2 = Pe;
tf = (zr < zf - lambda);
Pe2(tf) = 0;


%
h1 = M1b*Pe2;
h2 = M2b*Pe2;
close all
scalingFactor =  3.2563;
plot(zr,h1,'o',zr,h2,'x',z,h0/scalingFactor,'k');

%% calculate Pressure given h

h1b = h1;
h1b(nz) = (2^0.5);

Pe_check = abs(M1b(1:nz-1,1:nz-1))\h1b(1:nz-1);
Pe_check2 = M1b\h1b;
%Pe_check = M1invs*h1b(1:nz-1);
%Does this equal vector 1? No.

h2b = h2;
h2b(nz) = (2^0.5);

Pe_check2 = M2b\h2b;
close all
figure; hold on
plot(zr(1:nz-1),Pe_check,'o-r')
%plot(zr,Pe_check,'o-r')
%plot(zr,Pe_check2,'o-b')

%%
M1inv = inv(M1b);
M2inv = inv(M2b);

%inverse of symmetric matrices
M1s = M1b(1:nz-1,1:nz-1);
M1invs = inv(M1s);

M2s = M2b(1:nz-1,1:nz-1);

%check if matrix is positive definite
is_SPD(M1)

check = M1s*M1invs;


M2invs = inv(M2b(1:nz-1,1:nz-1));


%% approximate Pe using gmres

Kc = 0.1;
% h3 = h1b(1:nz-1);
% h3b = vertcat(h3,(Kc*2^0.5));
h1b(nz) = Kc*2^0.5;
tol = 1e-2;
maxit = 999;
%Pe_gm = gmres(M1s,h2b(1:nz-1),[],tol,maxit); 
L = ilu(M1b); %eye(size(M1b)
Pe_gm = gmres(M1b,h1b,[],tol,maxit,L); 

close all
plot(zr,Pe_gm,'o')

%%
close all
figure(1); pc = pcolor(M1b)
colorbar
pc.LineStyle='none'


%%

B = M1.'; %Tranpose of M1
X = M1-B;
%%

i = nz;
j = 2;

klogterm(zr(i),zr(j),zf)
klogterm(zr(j),zr(i),zf)

diff = zeros(nz-1,1);
i = 1;
for k = 1:nz-1
    diff(k) = M1b(i,k+1) - M1b(i,k);
end

%%
i = 1;
K = zeros(nz,1);
%K = log(abs( ((zf - zr(i))^0.5 + (zf - zr).^0.5)/ ( (zf - zr(i))^0.5 - (zf - zr).^0.5))) ; 

for j = 1:nz
    z0 = zr(j);
    z1 = zr(j) - dz/2;
    z2 = zr(j) + dz/2;
    if i == j
        K(j) = klogterm(zr(i),z1,zf) + klogterm(zr(i),z2,zf); 
    else
        K(j) = 2*klogterm(zr(i),z0,zf) + klogterm(zr(i),z1,zf) + klogterm(zr(i),z2,zf); 
    end
end


close all
figure; hold on
plot(zr,K)



function kij = klogterm(i,j,zf)
    kij = log( abs( ((zf - i)^0.5 + (zf - j)^0.5)/((zf - i)^0.5 - (zf - j)^0.5) ) ); 
end










