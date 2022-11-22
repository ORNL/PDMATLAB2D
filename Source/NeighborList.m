
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
% The function NeighborList generates a neighbor list and computes 
% per-bond quantities 
% ========================================================================

% Input
% -----
% Nx        : number of nodes in the x-direction
% Ny        : number of nodes in the y-direction
% xx        : x coordinates of all nodes in the grid (1D array of length Nx x Ny)
% yy        : y coordinates of all nodes in the grid (1D array of length Nx x Ny)
% xx1       : x coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% yy1       : y coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% M         : array that maps each node to its cell vertices
%             (2D array of size (Nx x Ny)-by-(number of vertices per cell))   
% del       : horizon
% dx        : mesh spacing in the x-direction
% dy        : mesh spacing in the y-direction
% VV        : area of each cell corresponding to each node (1D array of length Nx x Ny)
% omega     : influence function order indicator
%             0  : constant = 1
%             0.5: piecewise constant
%             1  : piecewise linear
%             3  : piecewise cubic
%             5  : piecewise quintic
%             7  : piecewise septic
% AlgName   : name of algorithm (in string format) for computation of neighbor areas
%             'FA'    : FA algorithm
%             'PA-AC' : PA-AC algorithm
%             'IPA-AC': IPA-AC algorithm
% flag_RDUG : if == 1, then the grid is uniform over a rectangular domain
%             (dx = dy); the flag name RDUG stands for "Rectangular Domain
%             Uniform Grid"
%             if == 0, then the grid is a general grid

% Output
% ------
% u_NA      : array of neighbor numbers for all nodes (2D array for most cases)
% IF_NA     : array of influence function values of neighbor bonds for all nodes (2D array for most cases)
% V_NA      : array of neighbor areas for all nodes (2D array for most cases)
% r_hat_NA  : array of reference lengths of neighbor bonds for all nodes (2D array for most cases)
% x_hat_NA  : array of x-coordinates of quadrature points for all nodes (2D array for most cases)
% y_hat_NA  : array of y-coordinates of quadrature points for all nodes (2D array for most cases)

% Discussion:
% ----------
% The following quantities are pre-computed per bond (for bond i-k with source node i):
% - Influence function (IFk)
% - neighbor area (Vk)
% - reference bond length (rk_hat)
% - coordinates of quadrature point (xk_hat and yk_hat)

% For the case of a uniform grid over a rectangular domain (flag_RDUG = 1), 
% the function employs an efficient algorithm to loop over the cells that
% overlap the neighborhood of a given source node. This algorithm was presented
% in the reference below (see Algorithm PA-AC Family Interaction in page 190).

% Reference: 
% ---------
% P. Seleson, Improved one-point quadrature algorithms for two-dimensional
% peridynamic models based on analytical calculations, Computer Methods in 
% Applied Mechanics and Engineering 282 (2014): 184–217.

function [u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG)

    % Tolerance
    tol = 1E-15;

    if flag_RDUG == 1

        % ----------------------------------------------------------------
        %                    Check validity of inputs
        % ----------------------------------------------------------------

        % Number of nodes
        Nnodes = Nx*Ny;

        % Check consistency of xx and yy array dimensions with number of nodes
        if length(xx) ~= Nnodes || length(yy) ~= Nnodes
            error('Nx and Ny are inconsistent with the dimensions of xx or yy arrays.')
        end

        % Number of vertices
        Nvertices = (Nx+1)*(Ny+1);

        % Check consistency of xx1 and yy1 array dimensions with number of vertices
        if length(xx1) ~= Nvertices || length(yy1) ~= Nvertices
            error('Nx and Ny are inconsistent with the dimensions of xx1 or yy1 arrays.')
        end

        % Size of M
        [s1,s2] = size(M);

        % Check consistency of M array dimensions with both number of nodes and rectangular grid
        if s1 ~= Nnodes || s2 ~= 4
            error('Nx and Ny are inconsistent with the dimensions of M, or M is inconsistent with a rectangular grid.')
        end

        % Check grid regularity
        xx1r = reshape(xx1,Nx+1,Ny+1);
        yy1r = reshape(yy1,Nx+1,Ny+1);

        maxdx = max(max(diff(xx1r,1,1)));
        mindx = min(min(diff(xx1r,1,1)));

        maxdy = abs(max(max(diff(yy1r,1,2))));
        mindy = abs(min(min(diff(yy1r,1,2))));

        if abs( maxdx - mindx ) > tol || abs( maxdy - mindy ) > tol
            error('Grid is not regular.')
        end

        % Check consistency of dx and dy with grid spacings
        if abs( maxdx - dx ) > tol || abs( maxdy - dy ) > tol
            error('dx or dy are inconsistent with the grid spacing in the x- or y-directions, respectively.')
        end

        % Check if the grid is uniform
        if abs(dx - dy) < tol 

            % Define uniform grid spacing
            h = dx;

            % ------------------------------------------------------------
            %            Check consistency between xx and xx1 
            %                and between yy and yy1 arrays
            % ------------------------------------------------------------

            % Loop through all nodes
            for n = 1:Nnodes

                % Define cell limits for node n
                minx = xx(n) - 0.5*h;
                maxx = xx(n) + 0.5*h;
                                
                miny = yy(n) - 0.5*h;
                maxy = yy(n) + 0.5*h;

                % Check if cell vertices are consistent with cell limits
                if max(abs(xx1(M(n,:)) - [minx maxx maxx minx]')) > tol || max(abs(yy1(M(n,:)) - [maxy maxy miny miny]')) > tol
                    error('xx or yy are inconsistent with xx1 or yy1, respectively.')
                end

            end

            % ------------------------------------------------------------
            %              Neighbor list for uniform grid
            %                 over rectangular domain
            % ------------------------------------------------------------

            % Initialization of variables
            % ---------------------------
            u_NA     = zeros(Nnodes,1); % Array of neighbor numbers for all nodes
            IF_NA    = zeros(Nnodes,1); % Array of influence function values of neighbor bonds for all nodes
            V_NA     = zeros(Nnodes,1); % Array of neighbor areas for all nodes
            r_hat_NA = zeros(Nnodes,1); % Array of reference lengths of neighbor bonds for all nodes
            x_hat_NA = zeros(Nnodes,1); % Array of x-coordinates of quadrature points for all nodes
            y_hat_NA = zeros(Nnodes,1); % Array of y-coordinates of quadrature points for all nodes

            % Compute maximum number of one-sided neighbor cells 
            % overlapping with the neighborhood of a source node
            Nn = floor(del/h + 0.5 - tol);

            % Loop through all nodes: (i,j)
            for i = 1:Ny    % Loop top to bottom
                for j = 1:Nx    % Loop left to right

                    % Convert (i,j) into node number
                    ui = (i-1)*Nx + j;

                    % Initializations
                    % 
                    % Note: We overestimate the number of neighbors by considering all cells within a square region 
                    %       containing the neighborhood of a source node. 
                    %       The "-1" in the initializations below subtracts the self node.
                    uVec     = zeros(1,((2*Nn+1)^2)-1); % Array of neighbor numbers for source node ui
                    IFVec    = zeros(1,((2*Nn+1)^2)-1); % Array of influence function values of neighbor bonds for source node ui
                    VVec     = zeros(1,((2*Nn+1)^2)-1); % Array of neighbor areas for source node ui
                    r_hatVec = zeros(1,((2*Nn+1)^2)-1); % Array of reference lengths of neighbor bonds for source node ui
                    x_hatVec = zeros(1,((2*Nn+1)^2)-1); % Array of x-coordinates of quadrature points for source node ui
                    y_hatVec = zeros(1,((2*Nn+1)^2)-1); % Array of y-coordinates of quadrature points for source node ui

                    % Counter initialization
                    counter = 0;

                    % Loop through potential neighbor nodes: (k,l) to check if their cells overlap with the neighborhood of node ui
                    for k = i - Nn:i + Nn % Loop top to bottom
                        
                        % Check if index k is within the range of possible values for the y-direction
                        if k < 1 || k > Ny

                        else

                            % For a given y location determined by k, calculate the maximum number of 
                            % one-sided neighbor cells in the horizontal direction overlapping with 
                            % the neighborhood of source node ui
                            if k == i
                                Nl = Nn;
                            else
                                % Find node above (for k < i) or below (for k > i) node i having indices (k,j)
                                uk1 = (k-1)*Nx + j;

                                xi_2 = abs(yy(uk1)-yy(ui));
                                Nl = floor( sqrt(del^2 - (xi_2 - 0.5*h)^2) / h  + 0.5 - tol);
                            end

                            for l = j - Nl:j + Nl % Loop left to right

                                % Check if index l is within the range of possible values for the x-direction
                                if l < 1 || l > Nx

                                else

                                    % Omit node (k,l) if it is the same as node (i,j)
                                    if i == k && j == l

                                    else
                                        % Convert (k,l) into node number
                                        uk = (k-1)*Nx + l;

                                        % Coordinates of four vertices of cell uk
                                        xx1k = xx1(M(uk,:));
                                        yy1k = yy1(M(uk,:));

                                        % Compute neighbor area, bond length, and coordinates of quadrature point
                                        [Vk,rk_hat,xk_hat,yk_hat] = NeighborAreaBondLengthCoord(xx(ui),xx(uk),yy(ui),yy(uk),dx,dy,xx1k,yy1k,del,VV(uk),AlgName);

                                        % Only include neighbor if neighbor area > 0
                                        if Vk > 0

                                            % Evaluate influence function
                                            [IFk] = InfluenceFunction(omega,rk_hat,del);

                                            % Update counter
                                            counter = counter + 1;

                                            % Node uk number
                                            uVec(counter) = uk;
                                            % Influence function value for bond ui-uk
                                            IFVec(counter) = IFk;
                                            % Neighbor area for bond ui-uk
                                            VVec(counter) = Vk;
                                            % Bond ui-uk reference length
                                            r_hatVec(counter) = rk_hat;
                                            % x-coordinate of quadrature point associated to bond ui-uk
                                            x_hatVec(counter) = xk_hat;
                                            % y-coordinate of quadrature point associated to bond ui-uk
                                            y_hatVec(counter) = yk_hat;

                                        end
                                    end
                                end
                            end
                        end
                    end

                    % Store computed per-bond variables for each node
                    u_NA(ui,1:counter)     = uVec(1:counter);
                    IF_NA(ui,1:counter)    = IFVec(1:counter);
                    V_NA(ui,1:counter)     = VVec(1:counter);
                    r_hat_NA(ui,1:counter) = r_hatVec(1:counter);
                    x_hat_NA(ui,1:counter) = x_hatVec(1:counter);
                    y_hat_NA(ui,1:counter) = y_hatVec(1:counter);

                end
            end

        else

            error('Grid is not uniform.')

        end

    elseif flag_RDUG == 0

        % ----------------------------------------------------------------
        %                    Check validity of inputs
        % ----------------------------------------------------------------

        % Check if xx and yy arrays have the same length
        if length(xx) ~= length(yy)
            error('Arrays xx and yy have different lengths.')
        end

        % Check if VV array has the same length as xx and yy arrays 
        if length(VV) ~= length(xx)
            error('Array VV has different length than xx and yy arrays.')
        end

        % ---------------------------------------------------------------
        %               Neighbor list for general grids
        % ---------------------------------------------------------------

        % Number of nodes
        Nnodes = length(xx);

        % Initialization of neighbor arrays
        u_NA     = zeros(Nnodes,1); % Array of neighbor numbers for all nodes
        IF_NA    = zeros(Nnodes,1); % Array of influence function values of neighbor bonds for all nodes
        V_NA     = zeros(Nnodes,1); % Array of neighbor areas for all nodes
        r_hat_NA = zeros(Nnodes,1); % Array of reference lengths of neighbor bonds for all nodes
        x_hat_NA = zeros(Nnodes,1); % Array of x-coordinates of quadrature points for all nodes
        y_hat_NA = zeros(Nnodes,1); % Array of y-coordinates of quadrature points for all nodes

        % Loop through all nodes ui:
        for ui = 1:Nnodes

            % Initializations
            uVec     = []; % Array of neighbor numbers for source node ui
            IFVec    = []; % Array of influence function values of neighbor bonds for source node ui
            VVec     = []; % Array of neighbor areas for source node ui
            r_hatVec = []; % Array of reference lengths of neighbor bonds for source node ui
            x_hatVec = []; % Array of x-coordinates of quadrature points for source node ui
            y_hatVec = []; % Array of y-coordinates of quadrature points for source node ui
            counter  = 0;

            % Loop through all nodes uk:
            for uk = 1:Nnodes

                % Omit node uk if it is the same as node ui
                if ui == uk

                else

                    % Check if M is empty
                    if isempty(M)
                        xx1k = [];
                        yy1k = [];
                    else
                        % Coordinates of vertices of cell uk
                        xx1k = xx1(M(uk,:));
                        yy1k = yy1(M(uk,:));
                    end

                    % Compute neighbor area, bond length, and coordinates of quadrature point
                    [Vk,rk_hat,xk_hat,yk_hat] = NeighborAreaBondLengthCoord(xx(ui),xx(uk),yy(ui),yy(uk),dx,dy,xx1k,yy1k,del,VV(uk),AlgName);

                    % Only include neighbor if neighbor area > 0
                    if Vk > 0

                        % Evaluate influence function
                        [IFk] = InfluenceFunction(omega,rk_hat,del);

                        % Update counter
                        counter = counter + 1;

                        % Node uk number
                        uVec(counter) = uk;
                        % Influence function value for bond ui-uk
                        IFVec(counter) = IFk;
                        % Neighbor area for bond ui-uk
                        VVec(counter) = Vk;
                        % Bond ui-uk reference length
                        r_hatVec(counter) = rk_hat;
                        % x-coordinate of quadrature point associated to bond ui-uk
                        x_hatVec(counter) = xk_hat;
                        % y-coordinate of quadrature point associated to bond ui-uk
                        y_hatVec(counter) = yk_hat;

                    end
                end
            end

            % Store computed per-bond variables for each node
            u_NA(ui,1:counter)     = uVec(1:counter);
            IF_NA(ui,1:counter)    = IFVec(1:counter);
            V_NA(ui,1:counter)     = VVec(1:counter);
            r_hat_NA(ui,1:counter) = r_hatVec(1:counter);
            x_hat_NA(ui,1:counter) = x_hatVec(1:counter);
            y_hat_NA(ui,1:counter) = y_hatVec(1:counter);

        end

    else

        error('flag_RDUG should be 0 or 1.')

    end
    
end
