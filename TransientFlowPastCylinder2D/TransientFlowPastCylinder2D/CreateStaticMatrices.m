function [Me_uu Ke_uu Ke_up Ke_pp Ke_pu fe_m fe_u fe_p] = CreateStaticMatrices(xElem,mu,tau,u0Elem,Ordr,dim)

    pOrdr = max(2*Ordr,2*(Ordr-1));

    %%% Gauss Quadrature Integration
    [xi eta w] = gaussQuad2DTriangular(pOrdr);

    %%% initialize matrices
    Me_uu = zeros(3*Ordr*dim,3*Ordr*dim);
    Ke_uu = zeros(3*Ordr*dim,3*Ordr*dim);
    Ke_up = zeros(3*Ordr*dim,3*Ordr);
    Ke_pp = zeros(3*Ordr,3*Ordr);
    Ke_pu = zeros(3*Ordr,3*Ordr*dim);
    fe_m  = zeros(3*Ordr*dim,1);
    fe_u  = zeros(3*Ordr*dim,1);
    fe_p  = zeros(3*Ordr,1);

    for gp=1:size(xi,2)

        J = Jacobian2DTriangular(xElem,xi(gp),eta(gp),Ordr);
        detJ = det(J);
        invJ = inv(J);

        N = basisFunc2DTriangular(xi(gp),eta(gp),Ordr);
        G = invJ * gradN2DTriangular(xi(gp),eta(gp),Ordr);

        [Nu Gu divN] = VeloctyShapeFunctions(N,G,dim);

        x = N*xElem;
        u0 = N*u0Elem;
        % g = [Source(x) Source(x)];
        g = [0 0];

        W  = N;
        Gw = G;

        Wu  = Nu;
        Guw = Gu;

        wDetJ = w(gp) * detJ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Element Mass Matrix %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Muu = Int( Nu' * Nu ) dx
        %     = Int( Nu' * Nu * det(J)) dxi
        Me_uu = Me_uu + ( Nu' * Nu ) * wDetJ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Element Stiffness Matrix %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Kuu = Int( dWu/dx' * k * dNu/dx ) dx
        %    = Int( Guw' * k * Gw * det(J)) dxi
        for j=1:size(Gu,3)
            for k=1:size(Gu,3)
                Ke_uu(j,k) = Ke_uu(j,k) + ( mu * sum(sum((Guw(:,:,j) .* Gu(:,:,k)))) ) * wDetJ;
            end
        end

        % Kpp = Int( tau * dN/dx' * dN/dx ) dx
        %     = Int( tau * G' * G * det(J)) dxi
        Ke_pp = Ke_pp + ( tau * Gw' * G ) * wDetJ;


        % Kup = -Int( dNu/dx_i' * N ) dx
        %     = -Int( divN' * N * det(J)) dxi
        Ke_up = Ke_up - ( divN' * N ) * wDetJ;


        % Kpu = Int( W' * dNu/dx_i ) dx
        %     = Int( W' * divN * det(J)) dxi
        Ke_pu = Ke_pu + ( W' * divN ) * wDetJ;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Element Force Vector %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fm = Int( Wu' * u0 )dx
        for j=1:2:size(Nu,2)
            fe_m(j+0) = fe_m(j+0) + ( Wu(1,j+0) * u0(1) ) * wDetJ;
            fe_m(j+1) = fe_m(j+1) + ( Wu(2,j+1) * u0(2) ) * wDetJ;
        end

        % fu = Int( Wu' * f )dx
        for j=1:2:size(Nu,2)
            fe_u(j+0) = fe_u(j+0) + ( Wu(1,j+0) * g(1) ) * wDetJ;
            fe_u(j+1) = fe_u(j+1) + ( Wu(2,j+1) * g(2) ) * wDetJ;
        end

        % fp = Int( tau * Gw' * f )dx
        fe_p = fe_p + ( tau * ( Gw(1,:)*g(1) + Gw(2,:)*g(2) )' ) * wDetJ;

    end

end