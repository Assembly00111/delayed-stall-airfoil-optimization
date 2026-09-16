clc
clear
% close all

%% ------------------------------------------------------------------------
% Problem Setup
% -------------------------------------------------------------------------
dim = 2;
mu = 1/200;
iterMax = 10000;

%% ------------------------------------------------------------------------
% Geometry
% -------------------------------------------------------------------------
model = createpde(1);
model = domainDefinition(model);

figure(1)
pdegplot(model,'EdgeLabels','on','FaceLabels','on');

%% ------------------------------------------------------------------------
% Mesh
% -------------------------------------------------------------------------
Hmax = 2.0;
Hmin = 0.001;
ordr = "linear";

if ordr=="linear"
    Ordr = 1;
elseif ordr=="quadratic"
    Ordr = 2;
else
    error("Unsupported element order")
end

model = generateMesh(model,...
    'Hmax',Hmax,...
    'Hmin',Hmin,...
    'Hface',{[3],0.1,[2],0.2},...
    'Hedge',{[12,13,14,15],0.02},...
    'GeometricOrder',ordr);

figure(2)
pdemesh(model)

[p,e,t] = meshToPet(model);

Coords  = p';
elemCon = t(1:end-1,:)';

% elemCon1 = [ elemCon(:,[1,4,6])
%              elemCon(:,[4,2,5])
%              elemCon(:,[6,5,3])
%              elemCon(:,[5,6,4]) ];

nElems = size(elemCon,1);
nNodes = size(Coords,1);

%% ------------------------------------------------------------------------
% Boundary Conditions
% -------------------------------------------------------------------------
u1DirichletBoundary = [4,12,13,14,15];
u2DirichletBoundary = [1,3,4,12,13,14,15];
pDirichletBoundary  = [2];

u1NeumannBoundary = [];
u2NeumannBoundary = [];

%% Dirichlet BC
[DirichletNodes1,DirichletValues1] = ...
    CalculateDirichletNodes(model,Coords,u1DirichletBoundary,1);

[DirichletNodes2,DirichletValues2] = ...
    CalculateDirichletNodes(model,Coords,u2DirichletBoundary,2);

[DirichletNodes3,DirichletValues3] = ...
    CalculateDirichletNodes(model,Coords,pDirichletBoundary,3);

DirichletIndex = [ ...
    2*DirichletNodes1-1 ...
    2*DirichletNodes2 ...
    DirichletNodes3 + 2*nNodes];

DirichletValue = [ ...
    DirichletValues1 ...
    DirichletValues2 ...
    DirichletValues3];


%% Neumann BC
[NeumannElements,~,NeumannElement,NeumannElementMap] = ...
    CalculateNeumannElements(model,Coords,elemCon,u1NeumannBoundary);

fu1 = ApplyNeumannBC(elemCon,Coords,mu,...
    NeumannElements,[],NeumannElement,NeumannElementMap);

[NeumannElements,~,NeumannElement,NeumannElementMap] = ...
    CalculateNeumannElements(model,Coords,elemCon,u2NeumannBoundary);

fu2 = ApplyNeumannBC(elemCon,Coords,mu,...
    NeumannElements,[],NeumannElement,NeumannElementMap);

fn = zeros(dim*nNodes,1);
fn(1:2:end) = fu1;
fn(2:2:end) = fu2;

%% Wall Boundary (computed once)
wallBoundary = [12,13,14,15];
[wallElements,~,wallElement,wallElementMap] = ...
    CalculateNeumannElements(model,Coords,elemCon,wallBoundary);

wallNodeList = findNodes(model,"region","Edge",wallBoundary);


%% ------------------------------------------------------------------------
% Global Time Stepping
% -------------------------------------------------------------------------
nTime = 2000;
dt = 0.05;

%% Initial Guess
u = zeros(nNodes,2);
p = zeros(nNodes,1);
x = zeros(3*nNodes,1);
u(:,1) = 1;

x(1:2:2*nNodes) = u(:,1);
x(2:2:2*nNodes) = u(:,2);
x(2*nNodes+1:3*nNodes) = p(:,1);


% %% Global Static Matrices
% % [Muu Kuu Kup Kpu Kpp fm fu fp] = StaticMatrix(Coords,elemCon,u,mu,Hmax,Ordr,dim);
% [Muu Kuu Kup Kpu Kpp fm fu fp] = StaticMatrixParallel(Coords,elemCon,u,mu,Hmax,Ordr,dim);

nLoad = 10;
for nt=1:nTime

    dtm = dt*(1.3^nLoad)
    nLoad = max(nLoad-1,0);

    if nt==1
        time(nt) = dtm;
    else
        time(nt) = time(nt-1) + dtm;
    end

    iterMax = min(10,1.3^(20-nLoad));

    u0 = u;
    p0 = p;
    
    %% Global Static Matrices
    [Muu Kuu Kup Kpu Kpp fm fu fp] = StaticMatrixParallel(Coords,elemCon,u0,mu,Hmax,Ordr,dim);

    %% Apply Neumann
    fu = fu + fn;
    
    Kuu_t = Kuu + Muu/dtm;
    % fu_t  = fu  + fm/dtm;
    fu_t  = fu  + Muu*x(1:2*nNodes,1)/dtm;
    
    %% Nonlinear Iteration
    [p, u, x] = CoreSteadyState(Coords,elemCon,u0,p0,mu,Kuu_t,Kup,Kpu,Kpp,fu_t,fp,DirichletIndex,DirichletValue,iterMax,Ordr,dim);
    
    %% Forces
    [C_D(nt),C_L(nt),C_p,C_s] = ...
        CalculateFluidForce(elemCon,Coords,mu,u,p,...
        wallElements,wallElement,wallElementMap);
    disp([C_D(nt) C_L(nt)])
    
    theta = linspace(0,360,length(wallNodeList)+1);
    theta(end) = [];
        
    figure(5)
    plot(theta,C_p(wallNodeList,1))
        
    figure(6)
    plot(theta,C_s(wallNodeList,1))

    figure(7)
    plot(time,C_D)
    
    figure(8)
    plot(time,C_L)

    drawnow

end


%% Final Plot
figure(10)
trisurf(elemCon,Coords(:,1),Coords(:,2),u(:,1),...
    'EdgeColor','none','facecolor','interp')
view(2)
axis equal
colormap jet
colorbar
