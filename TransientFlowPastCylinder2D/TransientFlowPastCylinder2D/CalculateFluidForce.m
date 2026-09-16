function [C_D C_L cp cs] = CalculateFluidForce(elemCon,Coords,mu,u,p,wallElements,ElementID,ElementMap)

    nNodes = size(Coords,1);

    cp = zeros(nNodes,1);   % local pressure coefficient
    cs = zeros(nNodes,1);   % local viscous coefficient

    fp = zeros(nNodes,2);   % local pressure force
    fs = zeros(nNodes,2);   % local viscous force

    % F_D = 0;    % total drag force
    % F_L = 0;    % total lift force

    Ordr = size(wallElements,2)-1;
    
    for i=1:size(wallElements,1)

        %--------------------------------------
        % Parent element data
        %--------------------------------------
        elemID    = ElementID(i);
        elemNodes = elemCon(elemID,:);
        xElem     = Coords(elemNodes,:);
        uElem     = u(elemNodes,:);
        pElem     = p(elemNodes);

        %--------------------------------------
        % Boundary edge data
        %--------------------------------------
        boundNodes = wallElements(i,:);
        xb = Coords(boundNodes,:);
        pb = p(boundNodes);
        elMap = ElementMap(i,:);

        %--------------------------------------
        % 1D Gauss quadrature on edge
        %--------------------------------------
        [r,w] = gaussQuad1D(Ordr);

        for gp = 1:length(r)
            
            %--------------------------------------
            % Map Gauss point to (xi,eta)
            %--------------------------------------
            xi  = 0;
            eta = 0;

            if elMap(1)==2
                xi = (1-r(gp))/2;
            elseif elMap(1)==3
                eta = (1-r(gp))/2;
            end

            if elMap(end)==2
                xi = (1+r(gp))/2;
            elseif elMap(end)==3
                eta = (1+r(gp))/2;
            end

            %--------------------------------------
            % Edge Jacobian
            %--------------------------------------
            J1 = Jacobian1D(xb(:,1),r(gp),Ordr);
            J2 = Jacobian1D(xb(:,2),r(gp),Ordr);

            detJb = sqrt(J1^2 + J2^2);
            wDetJb = w(gp) * detJb;

            tgt = [J1; J2] / detJb;
            nrm = [-tgt(2); tgt(1)];  % outward normal

            %--------------------------------------
            % 2D element mapping
            %--------------------------------------
            J = Jacobian2DTriangular(xElem,xi,eta,Ordr);
            invJ = inv(J);

            N  = basisFunc2DTriangular(xi,eta,Ordr);
            dN = gradN2DTriangular(xi,eta,Ordr);

            G = invJ * dN;   % physical gradients

            %--------------------------------------
            % Velocity gradient (2x2)
            %--------------------------------------
            grad_u = G * uElem;

            tau = mu * (grad_u + grad_u');

            %--------------------------------------
            % Pressure at GP (boundary interpolation)
            %--------------------------------------
            Nb = basisFunc1D(r(gp),Ordr);
            p_gp = Nb * pb;

            %--------------------------------------
            % normal coefficient of pressure
            %--------------------------------------
            cp(boundNodes,1) = cp(boundNodes,1) ...
                - Nb' * (-p_gp) * w(gp);

            %--------------------------------------
            % normal coefficient of shear (skin friction)
            %--------------------------------------
            cs(boundNodes,1) = cs(boundNodes,1) ...
                - Nb' * (tgt'*(tau*nrm)) * w(gp);


            %--------------------------------------
            % Assemble nodal forces
            %--------------------------------------
            fp(boundNodes,1) = fp(boundNodes,1) ...
                - Nb' * (-p_gp*nrm(1)) * wDetJb;

            fp(boundNodes,2) = fp(boundNodes,2) ...
                - Nb' * (-p_gp*nrm(2)) * wDetJb;

            fs(boundNodes,1) = fs(boundNodes,1) ...
                - Nb' * (tau(1,:)*nrm) * wDetJb;

            fs(boundNodes,2) = fs(boundNodes,2) ...
                - Nb' * (tau(2,:)*nrm) * wDetJb;

            %--------------------------------------
            % Total drag/lift
            %--------------------------------------
            % F_D = F_D + traction(1) * wDetJb;
            % F_L = F_L + traction(2) * wDetJb;

        end
    end

    F_D = sum(fp(:,1)) + sum(fs(:,1));
    F_L = sum(fp(:,2)) + sum(fs(:,2));

    % Drag/Lift coefficients (ρ=1, U=1, D=1 assumed)
    C_D = 2 * F_D;
    C_L = 2 * F_L;

end