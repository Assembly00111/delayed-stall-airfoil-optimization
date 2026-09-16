function [Muu Kuu Kup Kpu Kpp fm fu fp] = StaticMatrix(Coords,elemCon,u0,mu,Hmax,Ordr,dim)

    nNodes = size(Coords,1);
    nElems = size(elemCon,1);

    %% ------------------------------------------------------------------------
    % Global Matrices
    % -------------------------------------------------------------------------
    Muu = zeros(nNodes*dim,nNodes*dim);
    Kuu = zeros(nNodes*dim,nNodes*dim);
    Kup = zeros(nNodes*dim,nNodes);
    Kpu = zeros(nNodes,nNodes*dim);
    Kpp = zeros(nNodes,nNodes);
    
    fm  = zeros(nNodes*dim,1);
    fu  = zeros(nNodes*dim,1);
    fp  = zeros(nNodes,1);
    
    %% ------------------------------------------------------------------------
    % Assembly
    % -------------------------------------------------------------------------
    for elem = 1:nElems
    
        nodes = elemCon(elem,:);
        xElem = Coords(nodes,:);
        u0Elem = u0(nodes,:);
    
        % velocity DOF index
        uelemIndex = reshape([2*nodes-1;2*nodes],1,[]);
        pelemIndex = nodes;
    
        % element size
        xe = [xElem(1:3,:) ones(3,1)];
        area = det(xe);
    
        p1 = norm(xElem(2,:)-xElem(1,:));
        p2 = norm(xElem(3,:)-xElem(2,:));
        p3 = norm(xElem(1,:)-xElem(3,:));
        peri = p1+p2+p3;
    
        h = 4*area/peri;
    
        ue = mean(u0Elem)*mean(u0Elem)';
        tau = calculateTau(h,mu,ue);
    
        [Me_uu,Ke_uu,Ke_up,Ke_pp,Ke_pu,fe_m,fe_u,fe_p] = ...
            CreateStaticMatrices(xElem,mu,tau,u0Elem,Ordr,dim);
    
        Muu(uelemIndex,uelemIndex) = Muu(uelemIndex,uelemIndex) + Me_uu;
        Kuu(uelemIndex,uelemIndex) = Kuu(uelemIndex,uelemIndex) + Ke_uu;
        Kup(uelemIndex,pelemIndex) = Kup(uelemIndex,pelemIndex) + Ke_up;
        Kpp(pelemIndex,pelemIndex) = Kpp(pelemIndex,pelemIndex) + Ke_pp;
        Kpu(pelemIndex,uelemIndex) = Kpu(pelemIndex,uelemIndex) + Ke_pu;
    
        fm(uelemIndex) = fm(uelemIndex) + fe_m;
        fu(uelemIndex) = fu(uelemIndex) + fe_u;
        fp(pelemIndex) = fp(pelemIndex) + fe_p;
    end

end