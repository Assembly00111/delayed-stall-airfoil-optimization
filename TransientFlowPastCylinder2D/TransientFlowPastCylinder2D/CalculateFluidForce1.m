function [C_D, C_L, fp, fs] = CalculateFluidForce1( ...
        elemCon, Coords, mu, u, p, ...
        wallElements, ElementID, ElementMap)

% ------------------------------------------------------------
% Computes pressure and viscous traction forces on wall
%
% F = ∫ ( -p I + mu (grad u + grad u^T) ) n dΓ
% ------------------------------------------------------------

    nNodes = size(Coords,1);

    fp = zeros(nNodes,2);   % pressure contribution
    fs = zeros(nNodes,2);   % viscous contribution

    F_D = 0;
    F_L = 0;

    Ordr = size(wallElements,2) - 1;

    for i = 1:size(wallElements,1)

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
            % Traction
            %--------------------------------------
            traction = (-p_gp * eye(2) + tau) * nrm;

            wDetJb = w(gp) * detJb;

            %--------------------------------------
            % Assemble nodal forces
            %--------------------------------------
            fp(boundNodes,1) = fp(boundNodes,1) ...
                + Nb' * (-p_gp*nrm(1)) * wDetJb;

            fp(boundNodes,2) = fp(boundNodes,2) ...
                + Nb' * (-p_gp*nrm(2)) * wDetJb;

            fs(boundNodes,1) = fs(boundNodes,1) ...
                + Nb' * (tau(1,:)*nrm) * wDetJb;

            fs(boundNodes,2) = fs(boundNodes,2) ...
                + Nb' * (tau(2,:)*nrm) * wDetJb;

            %--------------------------------------
            % Total drag/lift
            %--------------------------------------
            F_D = F_D + traction(1) * wDetJb;
            F_L = F_L + traction(2) * wDetJb;

        end
    end

    % Drag/Lift coefficients (ρ=1, U=1, D=1 assumed)
    C_D = 2 * F_D;
    C_L = 2 * F_L;

end
