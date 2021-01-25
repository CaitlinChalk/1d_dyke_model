function J = fdJacobian(parameters,h,Pe)

%function to calculate analytical Jacobian of f (size n-3)
n = length(h);
dt = parameters.dt;
dz = parameters.dz;

alpha = (3*dt)/(16*(dz^2));

%Minv = inv(M1);

Jkk = zeros(n-3,1);

hk = h(3:n-2);
hkp = h(4:n-1);
hkm = h(2:n-3);
Pk = Pe(3:n-2);
Pkp = Pe(4:n-1);
Pkm = Pe(2:n-3);

Jkk(1) = 1 + 3*(dt/(2*dz))*(h(2)^2 * ((Pe(2) - Pe(1))/dz + 1) );
%Jkk(2:n-3) = 1 - 3*alpha*(((hk + hkp).^2)*(Pkp - Pk - dz) - (hk + hkm).^2*(Pk - Pkm - dz));

Jkp = -3*alpha*( ((hkp + hk).^2).*(Pkp - Pk - dz) );
Jkm = 3*alpha*( ((hk + hkm).^2).*(Pk - Pkm - dz) );
Jkk(2:n-3) = 1 + Jkp + Jkm;

%  Jkm = 3*alpha*( ((h(3:n-2) + h(2:n-3)).^2).*(Pe(3:n-2) - Pe(2:n-3) - dz) );
%  
%  Jkp = -3*alpha*( ((h(3:n-2) + h(2:n-3)).^2).*(Pe(3:n-2) - Pe(2:n-3) - dz) );
% % % % 
%  Jkk(2:n-3) = 1 + vertcat(Jkp,0) + vertcat(0,Jkm);
%  
 i = vertcat((1:n-3)',(2:n-3)',(2:n-4)');
 j = vertcat((1:n-3)',(1:n-4)',(3:n-3)');
 v = vertcat(Jkk, Jkm, Jkp(1:end-1));
 
% 
 J = sparse(i,j,v);

%Non-tridiagonal terms;

%  for i = 1:n-3
%      k = i+1;    
%      for j = 1:n-3
%          l = j+1;
%         jnw = [i-1, i, i+1]; %j not wanted
%         if ~ismember(j, jnw) 
%             J(i,j) = alpha*(h(k)*(Minv(k+1,l) - 2*Minv(k,l) + Minv(k-1,l)) + ...
%                      h(k+1)*(Minv(k+1,l) - Minv(k,l)) - h(k-1)*(Minv(k,l) - Minv(k-1,l)) );
%         end
%     end
% end

%boundary conditions?



        

