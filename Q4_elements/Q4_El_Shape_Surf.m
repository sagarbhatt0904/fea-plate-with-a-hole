function [NshapeS] = Q4_El_Shape_Surf(NES,xi)

NshapeS(1) = (1-xi)/2;
NshapeS(2) = (1+xi)/2;

%NshapeS = NshapeS';