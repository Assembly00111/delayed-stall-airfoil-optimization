function [Kuu fu] = ApplyDirichletBC(Kuu,fu,DirichletNodes,DirichletValues)
    
    for i=1:size(DirichletNodes,2)

        fu = fu - Kuu(:,DirichletNodes(i))*DirichletValues(i);
    
        Kuu(DirichletNodes(i),:) = 0;
        Kuu(:,DirichletNodes(i)) = 0;
        Kuu(DirichletNodes(i),DirichletNodes(i)) = 1;
    
        fu(DirichletNodes(i)) = DirichletValues(i);
    end

end