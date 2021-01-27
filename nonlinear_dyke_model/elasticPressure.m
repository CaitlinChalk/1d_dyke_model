function Pe = elasticPressure(parameters,h,Pe0,M1,M2)

%function to obtain the dyke pressure along depth from width
%input - h - array of width for 1 - n-1. h(n) = dzf = zf|t+1 - zf|t
%      - Pe0 - value of elastic Pressure at previous time step
%      - n size of h
%output - Pe - array of elastic pressure

%Step 1: define matrices M1 and M2 in h|t+1 = M1.*Pe|t+1 + dzf*M2*Pe|t
%Step 2: find Minv = the inverse of M1
%Step 3: calculate the elastic pressure Pe|t+1 = Minv*h|t+1 -
%        dzf*Minv.*V2, V2 = M2*Pe|t

%zf = z(n), or dzf + zf0?

Kc = parameters.Kc;

n = length(h);

%define V2
%V2 = M2(1:n-1,1:n-1)*Pe0(1:n-1);
V2 = M2*Pe0;
%% Elastic pressure

h2 = h;
h2(n) = Kc*(2^0.5);

%Pe = (Minv*h2' - h(n)*Minv*V2)';
%Pe = (M1\h(1:n-1) - h(n)*(M1\V2));

tol = 5e-2;
maxit = n-1;
Pe1 = gmresnomsg(M1,h2,[],tol,maxit); 
% if size(Pe) < n
%     Pe = vertcat(Pe,0);
% end
Pe2 = gmresnomsg(M1,V2,[],tol,maxit); 
Pe = Pe1 + h(n)*Pe2;


%Pe(tf) = 0;
%To do - calculate Pe(n)?

end

function x = gmresnomsg(varargin)
    [x,~] = gmres(varargin{:});
end
