function [fu] = ApplyNeumannBC(elemCon,Coords,k,NeumannElements,NeumannGradients,ElementID,ElementMap)

    Ordr = size(NeumannElements,2)-1;

    fu = zeros(size(Coords,1),1);
    
    for i=1:size(NeumannElements,1)

        elemNodes = elemCon(ElementID(i),:);
        boundNodes = NeumannElements(i,:);

        xb = Coords(NeumannElements(i,:),:);

        pOrdr = Ordr;

        %%% Gauss Quadrature Integration
        [r w] = gaussQuad1D(pOrdr);

        elMap = ElementMap(i,:);

        xi  = 0*r;
        eta = 0*r;

        if elMap(1)==2
            xi  = (1-r)/2;
        else if elMap(1)==3
            eta = (1-r)/2;
        end
        end

        if elMap(Ordr+1)==2
            xi  = (1+r)/2;
        else if elMap(Ordr+1)==3
            eta = (1+r)/2;
        end
        end
    
        for gp=1:size(xi,2)

            J1 = Jacobian1D(xb(:,1),r(gp),Ordr);
            J2 = Jacobian1D(xb(:,2),r(gp),Ordr);
            detJb = sqrt(J1*J1+J2*J2);

            Nb = basisFunc1D(r(gp),Ordr);
    
            N = basisFunc2DTriangular(xi(gp),eta(gp),Ordr);

            gradT = Nb*NeumannGradients(i,:)';
    
            W = N;
    
            wDetJb = w(gp) * detJb;
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Element Force Vector %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fe = Int( W' * k * dT/dx_i n_i )dx

            fu(elemNodes,1) = fu(elemNodes,1) + ( W' * k * gradT ) * wDetJb;
        end
    end

end