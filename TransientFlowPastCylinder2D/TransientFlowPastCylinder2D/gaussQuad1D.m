function [xi w] = gaussQuad1D(Ordr)
    if Ordr<2
        xi = [0];
        w  = [2];
    end

    if Ordr==2
        xi = [-1/sqrt(3)  1/sqrt(3)];
        w  = [1           1  ];
    end
end