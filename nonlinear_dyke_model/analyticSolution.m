function h0 = analyticSolution(x,lambda)

%produces analytic solution for fracture width which is kept open by a unit
%overpressure applied over length lambda at the tip

h0 = zeros(length(x));
fx = ( (x.^0.5 + lambda^0.5)./(x.^0.5 - lambda^0.5) );
ft = -2*x.*atanh( (lambda./x).^0.5 );
h0 = (1/pi)*(ft + lambda*log(abs(fx)) + 2*(lambda*x).^0.5);

%  L = lambda;
% % tol = 10e-01;
% % 
% % tf = ((x < lambda + tol) &  (x > lambda - tol));
% % x2 = x;
% % x2(tf) = 0;
% 
% fz = ( (x - lambda).^0.5 + lambda^0.5 )./( (x - lambda).^0.5 - lambda^0.5 );
% h0 = (1/pi)*( 2*(lambda*(x)).^0.5 - (x - lambda).*log(abs(fz)) );
% 
% 
% % for i = 1:length(x)
% %     h0(i) = constructh(x(i),L);
% % end


end

%  function hi = constructh(x2,L)
% 
% f1 = L^2 - x2^2 + 2*(L^0.5)*(x2^(3/2)) - 2*(L^3/2)*(x2^0.5);
% flog = log(abs((-x2^0.5 - L^0.5)/(-x2^0.5 + L^0.5))); 
% f2 = -4*L*x2 + 2*(L^0.5)*(x2^(3/2)) + 2*(L^(3/2))*(x2^0.5);
% fden = -2.*(x2^0.5)*(L^0.5) + L + x2;
% 
% hi = (f1*flog + f2)/fden;
% 
% end