function h = solveSystem(parameters, h0, Pe0, LS, V, xn, xn1,Kc)


[A,b] = dykewidthFDM4(parameters,h0,Pe0,LS,V,xn,xn1,Kc);

 dA = decomposition(A);
 h = dA\b;
% E = norm(A*h - b);

end

