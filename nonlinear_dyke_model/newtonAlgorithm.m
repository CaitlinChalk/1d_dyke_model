function [h,Pe,f,z] = newtonAlgorithm(parameters, h0, Pe0, z, zf0, dykewidthFDM, fdJacobian, tol, maxIts)

%Input - parameters - model parameters
%      - dykewidthFDM - function handle for nonlinear system
%      - h0 - initial state for h
%      - tol - convergence tolerance
%      - maxIt - max no. of iterations

%Output - h - updated dyke width
%       - f - function value
n = length(h0);
h = h0; %initial guess for h
Pe = Pe0; %initial guess for pressure

%calculate initial f
f = dykewidthFDM(parameters,h,h0,Pe,Pe0,z,zf0);


[M1,M2] = matrixM(z,zf0);
it = 0;
while norm(f) > tol && it < maxIts
    J = fdJacobian(parameters,h,Pe);
    delta1 = linearSolve(J,f); %linear solve
    %delta2 = solveJac(J,f);
    %check = delta1 - delta2;
    h(2:n-2) = h(2:n-2) + delta1;
    %h2 = h(2:n-1) - h(1:n-2);
    Pe = elasticPressure(parameters, h, Pe0, z, zf0, M1, M2); %update Pressure
    f = dykewidthFDM(parameters,h,h0,Pe,Pe0,z,zf0);
    h(n) = tipIncrement2(parameters,z,h,h0);    
    z = z + h(n);
    %fprintf('norm(f) at Newton step %d is %.4f \n', it, norm(f))
    it = it + 1;
end
