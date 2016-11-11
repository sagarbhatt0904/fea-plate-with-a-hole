function [Jac,detJ,Jhat] = Q4_El_Jacobian(NE,xi,eta,xy,DNshape)

Jac = zeros(2,2);

for i=1:NE
    Jac(1,1) = Jac(1,1) + DNshape(i,1)*xy(i,1);
    Jac(1,2) = Jac(1,2) + DNshape(i,1)*xy(i,2);
    Jac(2,1) = Jac(2,1) + DNshape(i,2)*xy(i,1);
    Jac(2,2) = Jac(2,2) + DNshape(i,2)*xy(i,2);
end

detJ = det(Jac);
Jhat = inv(Jac);
