function G = dNdxi1D(xi,Ordr)
    if Ordr == 1
        G = [-1/2   1/2];
    end

    if Ordr == 2
        G = [-1/2+xi    -2*xi   1/2+xi];
    end
end