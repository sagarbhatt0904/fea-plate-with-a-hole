function [XGS,WGS] = Q4_El_Gauss_Points_Surf(NGS)

if (NGS == 2)
    
    alf = sqrt(1/3);

    XGS(1,1) = -alf;
    XGS(2,1) = +alf;

    WGS(1) = 1;
    WGS(2) = 1;
    
else
    
    alf = sqrt(3/5);

    XGS(1,1) = -alf;
    XGS(2,1) = 0;
    XGS(3,1) = +alf;

    WGS(1) = 5/9;
    WGS(2) = 8/9;
    WGS(3) = 5/9;
    
end