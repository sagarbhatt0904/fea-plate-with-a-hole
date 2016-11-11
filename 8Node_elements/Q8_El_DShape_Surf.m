function [DNshapeS] = Q8_El_DShape_Surf(NES,xi)

DNshapeS(1) = xi - 1/2;
DNshapeS(2) = -2*xi;
DNshapeS(3) = -xi - 1/2;