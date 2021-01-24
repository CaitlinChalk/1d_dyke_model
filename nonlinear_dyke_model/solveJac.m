function delta = solveJac(J,f)

tol = 1e-2;
maxit = 100;
delta = gmresnomsg(J,-f,[],tol,maxit);

end

function x = gmresnomsg(varargin)
    [x,~] = gmres(varargin{:});
end