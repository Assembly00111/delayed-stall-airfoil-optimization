function [p u x] = CoreSteadyState(Coords,elemCon,u0,p0,mu,Kuu,Kup,Kpu,Kpp,fu,fp,DirichletIndex,DirichletValue,iterMax,Ordr,dim)

    % elemCon1 = [ elemCon(:,[1,4,6])
    %          elemCon(:,[4,2,5])
    %          elemCon(:,[6,5,3])
    %          elemCon(:,[5,6,4]) ];

    nNodes = size(Coords,1);

    u = u0;
    p = p0;

    x(1:2:2*nNodes,1)        = u0(:,1);
    x(2:2:2*nNodes,1)        = u0(:,2);
    x(2*nNodes+1:3*nNodes,1) = p0(:,1);

    %% ------------------------------------------------------------------------
    % Nonlinear Iteration
    % -------------------------------------------------------------------------
    for iter = 1:iterMax
    
        x0 = x;
    
        uc = u;
    
        [Kuuc,Kpuc,fuc,fpc] = ...
            DynamicMatrixParallel(Coords,elemCon,u,uc,mu,Ordr,dim);
    
        Kuu1 = Kuu + Kuuc;
        Kup1 = Kup;
        Kpu1 = Kpu + Kpuc;
        Kpp1 = Kpp;
    
        % A = [Kuu1 Kup1
        %      Kpu1 Kpp1];
        A = sparse(3*nNodes,3*nNodes);
        A(1:2*nNodes,1:2*nNodes) = Kuu1;
        A(1:2*nNodes,2*nNodes+1:end) = Kup1;
        A(2*nNodes+1:end,1:2*nNodes) = Kpu1;
        A(2*nNodes+1:end,2*nNodes+1:end) = Kpp1;
    
        b = [fu; fp];
    
        [A,b] = ApplyDirichletBCRapid(A,b,...
            DirichletIndex,DirichletValue);
    
        x = LinearSolver(A,b);
    
        % res(iter) = norm(x-x0)^2;
        res(iter) = max(abs(x-x0));
    
        %% Extract solution
        u0 = u;
        u  = [x(1:2:2*nNodes) x(2:2:2*nNodes)];
        p  = x(2*nNodes+1:3*nNodes);
    
        %% Plot velocity magnitude
        figure(3)
        trisurf(elemCon,Coords(:,1),Coords(:,2),u(:,1),...
            'EdgeColor','none')
        view(2)
        axis equal
        colormap jet
        colorbar
        drawnow
    
        figure(4)
        semilogy(res(1:iter))
    
        if res(iter) < 1e-10
            break
        end
    
    end

end