function [detJS] = Q4_El_Jacobian_Surf(NES,xi,xyS,DNshapeS)

dxdxi = 0;
dydxi = 0;

for i=1:NES
    dxdxi = dxdxi + DNshapeS(i)*xyS(i,1);
    dydxi = dydxi + DNshapeS(i)*xyS(i,2);    
end

detJS = sqrt( dxdxi*dxdxi + dydxi*dydxi );