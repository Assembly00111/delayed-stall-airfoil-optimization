function [Ke_uuc Ke_puc fe_uc fe_pc] = CreateConvectiveMatrices(xElem,uElem,mu,tau,Ordr,dim)

    pOrdr = max(2*Ordr);

    %%% Gauss Quadrature Integration
    [xi eta w] = gaussQuad2DTriangular(pOrdr);

    %%% initialize matrices
    Ke_uuc = zeros(3*Ordr*dim,3*Ordr*dim);
    Ke_puc = zeros(3*Ordr,3*Ordr*dim);
    fe_uc  = zeros(3*Ordr*dim,1);
    fe_pc  = zeros(3*Ordr,1);

    for gp=1:size(xi,2)

        J = Jacobian2DTriangular(xElem,xi(gp),eta(gp),Ordr);
        detJ = det(J);
        invJ = inv(J);

        N = basisFunc2DTriangular(xi(gp),eta(gp),Ordr);
        G = invJ * gradN2DTriangular(xi(gp),eta(gp),Ordr);

        [Nu Gu divN] = VeloctyShapeFunctions(N,G,dim);

        x = N*xElem;
        u = N*uElem;
        dudx = G*uElem;

        % g = [Source(x) Source(x)];
        g = [0 0];

        W  = N;
        Gw = G;

        Wu  = Nu;
        Guw = Gu;

        wDetJ = w(gp) * detJ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Element Stiffness Matrix %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uGu(:,:) = ( u(1)*Gu(1,:,:) + u(2)*Gu(2,:,:) );

        % Kuuc = Int( Wu' * u_j * dNu/dx_j ) dx
        %      = Int( Wu' * u_j * Gw_j * det(J)) dxi
        Ke_uuc = Ke_uuc + ( Wu' * uGu ) * wDetJ;

        % Kpuc = Int( tau * dNdx * u_j * dNu/dx_j ) dx
        %      = Int( tau * G' * u_j * Gw_j * det(J)) dxi
        Ke_puc = Ke_puc + ( tau * Gw' * uGu ) * wDetJ;


        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%% Element Force Vector %%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ududx = [( u(1)*dudx(1,1)+u(2)*dudx(2,1) )
                 ( u(1)*dudx(1,2)+u(2)*dudx(2,2) )];

        % fuc = Int( W' * u_j * dNudx_j )dx
        for j=1:2:size(Nu,2)
            fe_uc(j+0) = fe_uc(j+0) - ( Wu(1,j+0) * ududx(1) ) * wDetJ;
            fe_uc(j+1) = fe_uc(j+1) - ( Wu(2,j+1) * ududx(2) ) * wDetJ;
        end

        % fpc = Int( W' * u_j * dNudx_j )ds
        fe_pc = fe_pc - ( tau * ( Gw(1,:)*ududx(1) + Gw(2,:)*ududx(2) )' ) * wDetJ;

    end

end