
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
% The function BondBreaking breaks critically stretched bonds, unless they 
% are connected to a point in a no-fail zone
% ========================================================================

% Input
% -----
% xx          : x coordinates of all nodes in the grid 
% yy          : y coordinates of all nodes in the grid 
% v           : displacement of each node in the x-direction
% w           : displacement of each node in the y-direction
% so          : critical stretch
% u_NA        : array of neighbor numbers for all nodes
% r_hat_NA    : array of reference lengths of neighbor bonds for all nodes 
% x_hat_NA    : array of x-coordinates of quadrature points for all nodes
% y_hat_NA    : array of y-coordinates of quadrature points for all nodes 
% mask_nofail : Boolean array with a value of 1 for cells in a no-fail zone
%               and 0 for the other cells (recall each cell has a node 
%               associated to it which shares the same number)

% Output
% ------
% u_NA        : updated array of neighbor numbers for all nodes (after possible bond breaking)

% Discussion:
% ----------
% The no-fail zone is a region containing points for which bonds connected to
% them are prevented to fail.

function [u_NA] = BondBreaking(xx,yy,v,w,so,u_NA,r_hat_NA,x_hat_NA,y_hat_NA,mask_nofail)

    % Tolerance
    tol = 1E-15;

    % Number of nodes
    Nnodes = length(xx);

    % Find maximum number of neighbors any node can have
    zmax = length(u_NA(1,:));

    % Loop through all nodes ui
    for ui = 1:Nnodes

        % Check if node ui is in the no-fail zone
        if mask_nofail(ui) == 1 
            continue
        end

        % Get x-coordinate of node ui
        xi = xx(ui);
        % Get y-coordinate of node ui
        yi = yy(ui);

        % Get x-component of displacement of node ui
        vi = v(ui);
        % Get y-component of displacement of node ui
        wi = w(ui);

        % Loop over neighbors of node ui
        for z = 1:zmax

            % Get neighbor cell number
            uk = u_NA(ui,z);

            % Only consider neighbor cells with a larger number
            % than the source node to compute the bond pairwise force
            % and pairwise potential once per bond (this implies uk > 0)
            if uk > ui

                % Check if cell uk is in the no-fail zone
                if mask_nofail(uk) == 1
                    continue
                end

                % Get x-coordinate of quadrature point
                xk_hat = x_hat_NA(ui,z);
                % Get y-coordinate of quadrature point
                yk_hat = y_hat_NA(ui,z);

                % Get x-component of displacement of node uk
                vk = v(uk);
                % Get y-component of displacement of node uk
                wk = w(uk);

                % Get bond ui-uk reference length
                Rk  = r_hat_NA(ui,z);

                % Compute deformed bond components and length
                rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                rk2 = rxk^2 + ryk^2;       % current bond length squared

                % Note: the critical stretch condition is:
                %       s := (sqrt(rk2) - Rk)/Rk >= so
                %       This is equivalent to the condition
                %             rk2 >= Rk^2 * (so+1)^2
                %       because s = sqrt(rk2)/Rk - 1

                if rk2  > ((so+1)*Rk)^2 - tol

                    % Remove cell uk from neighbor list of node ui
                    u_NA(ui,z) = 0;

                    % Read neighbor list of node uk
                    uVeck = u_NA(uk,:);

                    % Find where cell ui is in the neighbor list of node uk
                    mask = (uVeck==ui);

                    % Remove cell ui from neighbor list of node uk
                    u_NA(uk,mask) = 0;

                end
            end
        end
    end
end
