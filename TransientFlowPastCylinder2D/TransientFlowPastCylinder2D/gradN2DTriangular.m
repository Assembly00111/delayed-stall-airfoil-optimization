function G = gradN2DTriangular(xi,eta,Ordr)
    if Ordr == 1
        G = [-1     1       0
             -1     0       1];
    end

    if Ordr == 2
        tau = 1-xi-eta;
        G = [1-4*tau    4*xi-1  0           4*(tau-xi)  4*eta   -4*eta
             1-4*tau    0       4*eta-1     -4*xi       4*xi    4*(tau-eta)];
    end
end