function u = DirichletFunction(x,dof)

    if x(1)==-10 && dof==1
        u = 1;
    else
        u = 0;
    end

end