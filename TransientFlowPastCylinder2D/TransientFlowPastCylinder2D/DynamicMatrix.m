function [Kuuc Kpuc fuc fpc] = DynamicMatrix(Coords,elemCon,u0,uc,mu,Ordr,dim)

    nElems = size(elemCon,1);
    nNodes = size(Coords,1);
    Kuuc = zeros(nNodes*dim,nNodes*dim);
    Kpuc = zeros(nNodes,nNodes*dim);
    fuc  = zeros(nNodes*dim,1);
    fpc  = zeros(nNodes,1);

    for i=1:nElems
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Element Node Coordiantes %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xElem = Coords(elemCon(i,:),:);
        u0Elem = u0(elemCon(i,:),:);
        uElem = uc(elemCon(i,:),:);
    
        elemIndex = [];
        for j=1:size(elemCon,2)
            uelemIndex(dim*j-1) = dim*elemCon(i,j)-1;
            uelemIndex(dim*j-0) = dim*elemCon(i,j)-0;
            pelemIndex(j)       = elemCon(i,j);
        end

        xe = zeros(3,3)+1;
        xe(:,1) = xElem([1,2,3],1);
        xe(:,2) = xElem([1,2,3],2);

        area = det(xe);

        p1 = sqrt((xe(2,1)-xe(1,1))*(xe(2,1)-xe(1,1)) + (xe(2,2)-xe(1,2))*(xe(2,2)-xe(1,2)));
        p2 = sqrt((xe(3,1)-xe(2,1))*(xe(3,1)-xe(2,1)) + (xe(3,2)-xe(2,2))*(xe(3,2)-xe(2,2)));
        p3 = sqrt((xe(1,1)-xe(3,1))*(xe(1,1)-xe(3,1)) + (xe(1,2)-xe(3,2))*(xe(1,2)-xe(3,2)));

        peri = p1+p2+p3;

        h = 4*area/peri;

        ue = mean(u0Elem)*mean(u0Elem)';

        tau = calculateTau(h,mu,ue);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Element Matrix and Vector %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Ke_uuc Ke_puc fe_uc fe_pc] = CreateConvectiveMatrices(xElem,uElem,mu,tau,Ordr,dim);
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Addition to Global Matrix %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Kuuc(uelemIndex, uelemIndex) = Kuuc(uelemIndex, uelemIndex) + Ke_uuc;
        Kpuc(pelemIndex, uelemIndex) = Kpuc(pelemIndex, uelemIndex) + Ke_puc;
        fuc(uelemIndex)              = fuc(uelemIndex)              + fe_uc;
        fpc(pelemIndex)              = fpc(pelemIndex)              + fe_pc;
    
    end

end