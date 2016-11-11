function [NshapeS] = Q8_El_Shape_Surf(NES,xi)

NshapeS(1) = ((xi-0)*(xi-1))/((-1-0)*(-1-1));
NshapeS(2) = ((xi+1)*(xi-1))/((0+1)*(0-1));
NshapeS(3) = ((xi+1)*(xi-0))/((1+1)*(0-1));

%NshapeS = NshapeS';