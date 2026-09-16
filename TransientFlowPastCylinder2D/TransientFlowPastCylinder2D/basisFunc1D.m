function N = basisFunc1D(xi,Ordr)
    if Ordr == 1
        N = [(1-xi)/2   (1+xi)/2];
    end

    if Ordr == 2
        N = [(-xi/2)*(1-xi)     (1-xi)*(1+xi)   (xi/2)*(1+xi)];
    end
end