function [Kuu, fu] = ApplyDirichletBCRapid(Kuu, fu, DirichletNodes, DirichletValues)

    nodes = DirichletNodes(:);
    vals  = DirichletValues(:);

    % Adjust RHS
    fu = fu - Kuu(:,nodes) * vals;

    % Zero rows and columns
    Kuu(nodes,:) = 0;
    Kuu(:,nodes) = 0;

    % Set diagonal entries
    idx = sub2ind(size(Kuu),nodes,nodes);
    Kuu(idx) = 1;

    % Set RHS values
    fu(nodes) = vals;

end