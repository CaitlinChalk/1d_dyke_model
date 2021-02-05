function Minv = updateMinv(M1,h,Kc)

h2 = h;
n = length(h);
h2(n) = Kc*(2^0.5);

tol = 1e-2;
maxit = n-1;

res = 20;
Pe = gmresnomsg(M1,h2,res,tol,maxit);
%Pe = gmresnomsg(M1,h2,[],tol,maxit);
%Pe = elasticPressure(parameters, h, Pe0, M1,M2); %initial guess for pressure
Minv = Pe/h2;

end

function x = gmresnomsg(varargin)
    [x,~] = gmres(varargin{:});
end