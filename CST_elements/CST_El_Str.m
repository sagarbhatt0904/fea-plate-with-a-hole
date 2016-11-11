function [eps,str] = CST_El_Str(ipstrn,xy,u,h,Y,nu)

ndof = 6;

Abar = [ 1 xy(1,1) xy(1,2); 1 xy(2,1) xy(2,2); 1 xy(3,1) xy(3,2) ];
A = det(Abar)/2;

B = (1/A/2)*[ xy(2,2)-xy(3,2) 0 xy(3,2)-xy(1,2) 0 xy(1,2)-xy(2,2) 0;
              0 xy(3,1)-xy(2,1) 0 xy(1,1)-xy(3,1) 0 xy(2,1)-xy(1,1);
              xy(3,1)-xy(2,1) xy(2,2)-xy(3,2) xy(1,1)-xy(3,1) ...
              xy(3,2)-xy(1,2) xy(2,1)-xy(1,1) xy(1,2)-xy(2,2) ];

if (ipstrn == 1)
  c = Y*(1-nu)/(1-2*nu)/(1+nu);
  C = c*[ 1 nu/(1-nu) 0; nu/(1-nu) 1 0; 0 0 (1-2*nu)/(1-nu)/2 ];
else
  c = Y/(1-nu)/(1+nu);
  C = c*[ 1 nu 0; nu 1 0; 0 0 (1-nu)/2 ];
end

eps = B*u;
str = C*eps;





