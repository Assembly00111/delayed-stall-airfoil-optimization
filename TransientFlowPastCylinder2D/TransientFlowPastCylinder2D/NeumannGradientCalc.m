function ElemGrad = NeumannGradientCalc(xi)

    for i=1:size(xi,1)
        x = xi(i,1);
        y = xi(i,2);

        ElemGrad(i) = 0;

        if (x==1)
            ElemGrad(i) =  pi*cos(pi*x)*sin(pi*y);
        end
        if (x==0)
            ElemGrad(i) = -pi*cos(pi*x)*sin(pi*y);
        end
        if (y==1)
            ElemGrad(i) =  pi*sin(pi*x)*cos(pi*y);
        end
        if (y==0)
            ElemGrad(i) = -pi*sin(pi*x)*cos(pi*y);
        end

    end

end