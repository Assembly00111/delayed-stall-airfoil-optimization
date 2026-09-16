function [NeumannElements NeumannGradients ElementID ElementMap] = CalculateNeumannElements(model,Coords,elemCon,NeumannBoundary);

    Ordr = size(elemCon,2)/3;

    NeumannElements  = [];
    NeumannGradients = [];

    for i=1:size(NeumannBoundary,2)
        nodeList  = findNodes(model,"region","Edge",NeumannBoundary(i));
    
        ElemList = [];
        ElemGrad = [];

        k=1;
        for j=1:Ordr:size(nodeList,2)-1

            xi = Coords(nodeList(j:j+Ordr),:);

            ElemList(k,:) = nodeList(j:j+Ordr);
            ElemGrad(k,:) = NeumannGradientCalc(xi);
            k = k+1;
        end
    
        NeumannElements  = [NeumannElements  ;ElemList];
        NeumannGradients = [NeumannGradients ;ElemGrad];
    end

    ElementID = [];
    ElementMap = [];
    for i=1:size(NeumannElements,1)
        Elements1 = findElements(model,"attached",NeumannElements(i,1));
        Elements2 = findElements(model,"attached",NeumannElements(i,end));

        ElementID(i,1) = intersect(Elements1, Elements2);

        ElemNodes = elemCon(ElementID(i,1),:);

        idx1 = find(ElemNodes == NeumannElements(i,1));
        idx2 = find(ElemNodes == NeumannElements(i,2));

        if Ordr==1
            ElementMap(i,:) = [idx1 idx2];
        end
        if Ordr==2
            idx3 = find(ElemNodes == NeumannElements(i,3));
            ElementMap(i,:) = [idx1 idx2 idx3];
        end

    end

end