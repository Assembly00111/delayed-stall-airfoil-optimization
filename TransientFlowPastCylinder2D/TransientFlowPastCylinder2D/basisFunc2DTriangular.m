function N = basisFunc2DTriangular(xi,eta,Ordr)
    if Ordr == 1
        N = [(1-xi-eta)     xi      eta];
    end

    if Ordr == 2
        tau = 1-xi-eta;
        N = [tau*(2*tau-1)  xi*(2*xi-1)     eta*(2*eta-1)   4*tau*xi    4*xi*eta    4*tau*eta];
    end
end