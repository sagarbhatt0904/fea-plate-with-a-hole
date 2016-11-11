function [DNshape] = Q8_El_DShape(NE,xi,eta)


DNshape(:,1)=[- (xi/4 - 1/4)*(eta - 1) - ((eta - 1)*(eta + xi + 1))/4;
              ((eta - 1)*(eta - xi + 1))/4 - (xi/4 + 1/4)*(eta - 1);
              (xi/4 + 1/4)*(eta + 1) + ((eta + 1)*(eta + xi - 1))/4;
              (xi/4 - 1/4)*(eta + 1) + ((eta + 1)*(xi - eta + 1))/4;
                                                  xi*(eta - 1);
                                                1/2 - eta^2/2;
                                                 -xi*(eta + 1);
                                               eta^2/2 - 1/2]';
   

DNshape(:,2)=[- (xi/4 - 1/4)*(eta - 1) - (xi/4 - 1/4)*(eta + xi + 1);
             (xi/4 + 1/4)*(eta - xi + 1) + (xi/4 + 1/4)*(eta - 1);
             (xi/4 + 1/4)*(eta + 1) + (xi/4 + 1/4)*(eta + xi - 1);
             (xi/4 - 1/4)*(xi - eta + 1) - (xi/4 - 1/4)*(eta + 1);
                                               xi^2/2 - 1/2;
                                                -eta*(xi + 1);
                                               1/2 - xi^2/2;
                                                eta*(xi - 1)]';

