%=========================================================%
%                     DEMO_FE_MAN                         %
%=========================================================%


%---------------------------------------------------------%
% 1. Direct declaration of geometry                       %
% See section 3.1.1 of the tutorial                       %
%---------------------------------------------------------%
femesh('reset');
model = struct('Node',[],'Elt',[]);
%           NodeID  unused   x y z
model.Node=[  1      0 0 0    0 1 0;
              2      0 0 0    0 0 0;
              3      0 0 0    1 1 0;
              4      0 0 0    1 0 0;
              5      0 0 0    2 0 0;
              6      0 0 0    2 1 0;
              7      0 0 0    1 1 1]; % reference node

model.Elt=[
            % declaration of element group for longerons
                Inf     abs('beam1') 
            %node1  node2   MatID ProID nodeR, zeros to fill the matrix 
                1       3      1    1     7       0
                3       6      1    1     7       0
                2       4      1    1     7       0
                4       5      1    1     7       0
             % declaration of element group for diagonals
                Inf     abs('beam1')
                2       3      1    2     7       0
                4       6      1    2     7       0
             % declaration of element group for battens
                Inf     abs('beam1')
                3       4      1    3     7       0
                5       6      1    3     7       0 ];

feplot(model);
fecom('view2');
