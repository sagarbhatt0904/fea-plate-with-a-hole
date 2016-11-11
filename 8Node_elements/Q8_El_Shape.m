function [Nshape] = Q8_El_Shape(NE,xi,eta)


Nshape = [-1/4*(1-xi)*(1-eta)*(xi+eta+1);
          1/4*(1+xi)*(1-eta)*(xi-eta-1);
          1/4*(1+xi)*(1+eta)*(xi+eta-1);
         -1/4*(1-xi)*(1+eta)*(xi-eta+1);
                1/2*(1-xi^2)*(1-eta);
                1/2*(1-eta^2)*(1+xi);
                1/2*(1-xi^2)*(1+eta);
               1/2*(1-eta^2)*(1-xi)]'; 