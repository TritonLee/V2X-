function [x,fval] = goldmax(f,z1,z2,eps)
%%  nonlinear programming of Golden section research
 % designed to locate unimodal continuously functions 
    k = 0;
    g = (3-sqrt(5))/2;
    Z1 = z1+g*(z2-z1);
    Z2 = z2-g*(z2-z1);
    F1 = f(Z1);
    F2 = f(Z2);

    while abs(z2-z1)>=eps
        if  F1>F2
            z2 = Z2;
            Z2 = Z1;
            F2 = F1;
            Z1 = z1+g*(z2-z1);
            F1 = f(Z1);

        elseif F1<F2
            z1 = Z1;
            Z1 = Z2;
            F1 = F2;
            Z2 = z2-g*(z2-z1);
            F2 = f(Z2);
        else
            z1 = Z1;
            z2 = Z2;
            Z1 = z1+g*(z2-z1);
            Z2 = z2-g*(z2-z1);
            F1 = f(Z1);
            F2 = f(Z2);
        end
    end
    x = (z1+z2)/2;
    fval = f(x);


end



