function [XG,WG] = Q4_El_Gauss_Points(NG)

if (NG == 4)
    
    alf = sqrt(1/3);

    XG(1,1) = -alf;
    XG(2,1) = +alf;
    XG(3,1) = +alf; 
    XG(4,1) = -alf;

    XG(1,2) = -alf;
    XG(2,2) = -alf;
    XG(3,2) = +alf;
    XG(4,2) = +alf;

    for i=1:NG
        WG(i) = 1;
    end
    
else
    
    alf = sqrt(3/5);

    XG(1,1) = -alf;
    XG(2,1) = 0;
    XG(3,1) = +alf;
    XG(4,1) = -alf;
    XG(5,1) = 0;
    XG(6,1) = +alf;
    XG(7,1) = -alf;
    XG(8,1) = 0;
    XG(9,1) = +alf;
    
    XG(1,2) = -alf;
    XG(2,2) = -alf;
    XG(3,2) = -alf;
    XG(4,2) = 0;
    XG(5,2) = 0;
    XG(6,2) = 0;
    XG(7,2) = +alf;
    XG(8,2) = +alf;
    XG(9,2) = +alf;

    WG(1) = 25/81;
    WG(2) = 40/81;
    WG(3) = 25/81;
    WG(4) = 40/81;
    WG(5) = 64/81;
    WG(6) = 40/81;
    WG(7) = 25/81;
    WG(8) = 40/81;
    WG(9) = 25/81;
    
end
