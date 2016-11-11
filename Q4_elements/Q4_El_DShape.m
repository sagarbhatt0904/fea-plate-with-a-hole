function [DNshape] = Q4_El_DShape(NE,xi,eta)

DNshape(1,1) = -(1-eta)/4;
DNshape(2,1) = +(1-eta)/4;
DNshape(3,1) = +(1+eta)/4;
DNshape(4,1) = -(1+eta)/4;

DNshape(1,2) = -(1-xi)/4;
DNshape(2,2) = -(1+xi)/4;
DNshape(3,2) = +(1+xi)/4;
DNshape(4,2) = +(1-xi)/4;
