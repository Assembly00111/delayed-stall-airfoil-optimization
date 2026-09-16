function model = domainDefinition(model)

%%% Defining Domain Extends
R_1 = [ -10   -10
         25  -10
         25   10
        -10    10];

R_2 = [ -2   -2
         7   -3
         7    3
        -2    2];

R_3 = [  7   -3
        24.5   -3
        24.5    3
         7    3];

%%% Defining Cylinder Boundary
C_1 = [0      0             % centre of circle
       0.5    0   ];        % radius of circle

%%% generating model
R1 = [2,4,R_1(1,1),R_1(2,1),R_1(3,1),R_1(4,1),R_1(1,2),R_1(2,2),R_1(3,2),R_1(4,2)]';
R2 = [2,4,R_2(1,1),R_2(2,1),R_2(3,1),R_2(4,1),R_2(1,2),R_2(2,2),R_2(3,2),R_2(4,2)]';
R3 = [2,4,R_3(1,1),R_3(2,1),R_3(3,1),R_3(4,1),R_3(1,2),R_3(2,2),R_3(3,2),R_3(4,2)]';
C1 = [1,C_1(1,1),C_1(1,2),C_1(2,1)]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
gm = [R1,R2,R3,C1];
sf = 'R1-C1';  % boolean operation

ns = char('R1','R2','R3','C1');
ns = ns';
[dl bt] = decsg(gm,sf,ns);

% [g1,g2] = csgdel(dl,bt,[5]);
[g1,g2] = csgdel(dl,bt,[]);

geometryFromEdges(model,g1);