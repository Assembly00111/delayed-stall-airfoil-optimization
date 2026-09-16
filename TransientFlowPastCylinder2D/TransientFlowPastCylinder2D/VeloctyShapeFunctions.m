function [Nu Gu divN] = VeloctyShapeFunctions(N,G,dim)

    nNodesPerElem = size(N,2);

    % Nu = [  N(1)    0       N(2)    0       N(3)    0
    %         0       N(1)    0       N(2)    0       N(3)];
    % 
    % Gu(1,:,:) = [   G(1,1)  0       G(1,2)  0       G(1,3)  0
    %                 0       G(1,1)  0       G(1,2)  0       G(1,3)];
    % 
    % Gu(2,:,:) = [   G(2,1)  0       G(2,2)  0       G(2,3)  0
    %                 0       G(2,1)  0       G(2,2)  0       G(2,3)];
    % 
    % divN = [G(1,1)  G(2,1)  G(1,2)  G(2,2)  G(1,3)  G(2,3)];

    Nu = zeros(dim,dim*nNodesPerElem);
    Gu = zeros(dim,dim,dim*nNodesPerElem);
    divN = zeros(1,dim*nNodesPerElem);

    for i=1:size(N,2)
        for j=1:dim
            Nu(j,(i-1)*dim+j) = N(i);
            Gu(:,j,(i-1)*dim+j) = G(:,i);
        end
        divN(1,(i-1)*dim+1) = G(1,i);
        divN(1,(i-1)*dim+2) = G(2,i);
    end

end