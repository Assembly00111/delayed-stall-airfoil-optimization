function [DirichletNodes DirichletValues] = CalculateDirichletNodes(model,Coords,DirichletBoundary,dof)
    
    DirichletNodes = [];
    DirichletValues = [];
    
    for i=1:size(DirichletBoundary,2)
    
        nodeList  = findNodes(model,"region","Edge",DirichletBoundary(i));
        nodeValue = 0*nodeList;
    
        for j=1:size(nodeList,2)
            xi = Coords(nodeList(j),:);
            nodeValue(j) = DirichletFunction(xi,dof);
        end
    
        DirichletNodes  = [DirichletNodes  nodeList ];
        DirichletValues = [DirichletValues nodeValue];
    
    end

    [~, UniqueNodes] = unique(DirichletNodes, 'stable');

    DirichletNodes  = DirichletNodes(UniqueNodes);
    DirichletValues = DirichletValues(UniqueNodes);

    % size(DirichletNodes)

    if size(DirichletNodes,1)==0
        DirichletNodes  = [1];
        DirichletValues = [0];
    end

end