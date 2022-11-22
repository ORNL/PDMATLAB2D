
% ========================================================================
% Copyright (c) 2022 by Oak Ridge National Laboratory                      
% All rights reserved.                                                     
%                                                                           
% This file is part of PDMATLAB2D. PDMATLAB2D is distributed under a           
% BSD 3-clause license. For the licensing terms see the LICENSE file in    
% the top-level directory.                                                 
%                                                                          
% SPDX-License-Identifier: BSD-3-Clause                                    
% ========================================================================

% ========================================================================
% The function TestCriticalStretch checks the critical stretch expressions
% ========================================================================

% Discussion:
% ----------
% We compute a discrete approximation to the total energy per unit area 
% resulting from all the interactions between points along a horizontal line 
% to the left of a vertical line (representing a surface) and their neighboring
% points on the right side of the vertical line, under the assumption that 
% all the corresponding bonds are critically stretched. 
%
% The horizontal line is discretized with a set of segments of length h,
% while the right side of the vertical line is discretized with a uniform
% grid of cells of area h^2.
%
% We sum all pairwise potentials of bonds connecting the nodes centered at
% the line segments to their neighboring cells, weighted by the segment
% length times the neighboring area.
%
% The function computes values for two plane elasticity models:
% (i)  plane strain
% (ii) plane stress
% for different influence functions and neighbor areas algorithms and outputs
% the computed values along with the reference fracture energy value.

function TestCriticalStretch

    % Current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                    Material constants and model 
    % --------------------------------------------------------------------

    % Horizon
    del = 1;

    % Young's modulus
    E = 1;

    % Fracture energy
    Go = 1;

    % Model
    model = 'GPMB';

    % --------------------------------------------------------------------
    %                     Domain and discretization
    % --------------------------------------------------------------------
    
    % m-ratio
    m = 6;
    
    % Uniform grid spacing
    h = del/m;

    % Coordinates of intersection point between horizontal and vertical lines
    xI = 0;
    yI = 0;

    % --------------------------------------------------------------------
    %                          Create grids
    % -------------------------------------------------------------------- 

    % Grid spacings
    dx = h;
    dy = dx;
    
    % Left-side grid of nodes
    xxL  = xI-0.5*dx:-dx:xI-del+0.5*dx;
    dimL = length(xxL); 
    yyL  = yI*ones(dimL,1);
    
    % Right-side grid of nodes
    xR = xI+0.5*dx:dx:xI+del-0.5*dx;
    yR = yI-del:dy:yI+del;
    
    [xxR,yyR] = meshgrid(xR,yR);
    
    xxR  = xxR(:);
    yyR  = yyR(:);
    dimR = length(xxR);

    % --------------------------------------------------------------------
    %                 Compute total energy per unit area
    % -------------------------------------------------------------------- 

    % Run over plane elasticity models: plane strain and plane stress
    for PlanarModel = ["PlaneStrain" "PlaneStress"]

        if strcmp(PlanarModel,'PlaneStrain')
            fprintf('------------------------------------------------------- \n')
            fprintf('                    Plane strain \n')
            fprintf('------------------------------------------------------- \n')
        else
            fprintf('------------------------------------------------------- \n')
            fprintf('                    Plane stress \n')
            fprintf('------------------------------------------------------- \n')
        end

        % Print output labels
        fprintf('omega     AlgName       W         Go         \n')

        % Run over influence functions
        for omega = [0 0.5 1 3 5 7]

            % Run over neighbor areas algorithms
            for AlgName = ["FA" "PA-AC" "IPA-AC"]

                % PD material constants
                [c,so] = PDBondConstants(omega,del,E,Go,model,PlanarModel);

                % Initialize energy density
                W = 0;

                % Loop over nodes on horizontal line
                for nL = 1:dimL

                    % Get x-coordinate of node on horizontal line
                    xi = xxL(nL);
                    % Get y-coordinate of node on horizontal line
                    yi = yyL(nL);

                    % Loop over nodes on right side of vertical line
                    for nR = 1:dimR

                        % Get x-coordinate of node on right side of vertical line
                        xk = xxR(nR);
                        % Get y-coordinate of node on right side of vertical line
                        yk = yyR(nR);

                        % Coordinates of cell for node on right side of vertical line
                        minxk = xk - 0.5*h;
                        maxxk = xk + 0.5*h;
                        minyk = yk - 0.5*h;
                        maxyk = yk + 0.5*h;

                        xx1k = [minxk maxxk maxxk minxk];
                        yy1k = [maxyk maxyk minyk minyk];

                        % Area (in reference configuration) of cell for node on right side of vertical line
                        VVk = h^2;

                        % Neighbor area and reference bond length
                        [Vk,rk_hat,~,~] = NeighborAreaBondLengthCoord(xi,xk,yi,yk,dx,dy,xx1k,yy1k,del,VVk,AlgName);
                        Rk = rk_hat;

                        % Influence function
                        IFk = InfluenceFunction(omega,Rk,del);

                        % Compute bond pairwise potential under critical stretch so
                        wk = 0.5*c*IFk*(so^2)*Rk;

                        % Update energy density
                        W = W + wk*Vk*dx;

                    end

                end

                fprintf('  %3g      %6s     %5.3f     %5.3f \n',omega,AlgName,W,Go)

            end
        end
    end

    % Change directory
    cd(Directory)

end
