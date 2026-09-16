function [Kuuc,Kpuc,fuc,fpc] = DynamicMatrixParallel(Coords,elemCon,u0,uc,mu,Ordr,dim)

nElems = size(elemCon,1);
nNodes = size(Coords,1);

Kuuc = zeros(nNodes*dim,nNodes*dim);
Kpuc = zeros(nNodes,nNodes*dim);
fuc  = zeros(nNodes*dim,1);
fpc  = zeros(nNodes,1);

%% Storage for element contributions
Keuu_cell = cell(nElems,1);
Kepu_cell = cell(nElems,1);
feu_cell  = cell(nElems,1);
fep_cell  = cell(nElems,1);

uelemIndex_cell = cell(nElems,1);
pelemIndex_cell = cell(nElems,1);

%% ------------------------------------------------------------------------
% Parallel element computation
% -------------------------------------------------------------------------
parfor elem = 1:nElems

    nodes = elemCon(elem,:);
    xElem = Coords(nodes,:);
    u0Elem = u0(nodes,:);
    uElem  = uc(nodes,:);

    % DOF mapping
    uelemIndex = reshape([dim*nodes-1; dim*nodes],1,[]);
    pelemIndex = nodes;

    % Element size
    xe = [xElem(1:3,:) ones(3,1)];
    area = det(xe);

    p1 = norm(xElem(2,:)-xElem(1,:));
    p2 = norm(xElem(3,:)-xElem(2,:));
    p3 = norm(xElem(1,:)-xElem(3,:));
    peri = p1+p2+p3;

    h = 4*area/peri;

    ue = mean(u0Elem)*mean(u0Elem)';
    tau = calculateTau(h,mu,ue);

    % Element matrices
    [Ke_uuc,Ke_puc,fe_uc,fe_pc] = ...
        CreateConvectiveMatrices(xElem,uElem,mu,tau,Ordr,dim);

    % Store
    Keuu_cell{elem} = Ke_uuc;
    Kepu_cell{elem} = Ke_puc;
    feu_cell{elem}  = fe_uc;
    fep_cell{elem}  = fe_pc;

    uelemIndex_cell{elem} = uelemIndex;
    pelemIndex_cell{elem} = pelemIndex;
end

%% ------------------------------------------------------------------------
% Global assembly
% -------------------------------------------------------------------------
for elem = 1:nElems

    uelemIndex = uelemIndex_cell{elem};
    pelemIndex = pelemIndex_cell{elem};

    Kuuc(uelemIndex,uelemIndex) = Kuuc(uelemIndex,uelemIndex) + Keuu_cell{elem};
    Kpuc(pelemIndex,uelemIndex) = Kpuc(pelemIndex,uelemIndex) + Kepu_cell{elem};
    fuc(uelemIndex)             = fuc(uelemIndex)             + feu_cell{elem};
    fpc(pelemIndex)             = fpc(pelemIndex)             + fep_cell{elem};

end

end