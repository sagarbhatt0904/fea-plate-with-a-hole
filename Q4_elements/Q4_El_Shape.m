function [Nshape] = Q4_El_Shape(NE,xi,eta)

Nshape(1) = (1-xi)*(1-eta)/4;
Nshape(2) = (1+xi)*(1-eta)/4;
Nshape(3) = (1+xi)*(1+eta)/4;
Nshape(4) = (1-xi)*(1+eta)/4;
