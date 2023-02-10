
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
% The function PlotComplexShapes plots domains with shapes more complex 
% than a rectangular domain along with the grid and bonds
% ========================================================================

% Input
% -----
% shape          : domain shape
%                  'L-shape'       : L-shape domain
%                  'Circle'        : circular domain
%                  'SquareWithHole': square with a circular hole at its center
% flag_PlotBonds : if == 1, then the function plots bonds in addition to nodes
%                  if == 0, then the function does not plot bonds

function PlotComplexShapes(shape,flag_PlotBonds)

    % Identify current directory
    Directory = pwd;

    % Influence function order indicator
    omega = 0;

    % Algorithm for computation of neighbor areas
    AlgName = 'FA';

    % Flag for Rectangular Domain Uniform Grid (RDUG)
    flag_RDUG = 0;

    % --------------------------------------------------------------------
    %                        L-shape domain
    % --------------------------------------------------------------------

    if strcmp(shape,'L-shape')

        % Change directory
        cd ../Source/

        % ----------------------------------------------------------------
        %                      Generate subgrids
        % ----------------------------------------------------------------

        % Generate 1st uniform grid (grid a): bottom part of L-shape - rectangular domain (0,2)x(0,1)
        [xxa,yya,~,~,dxa,dya,VVa,xx1a,yy1a,Ma] = GridGenerator(0,2,0,1,10,10,0);

        % Generate 2nd uniform grid (grid b): top part of L-shape - square domain (0,1)x(1,2)
        [xxb,yyb,~,~,dxb,dyb,VVb,xx1b,yy1b,Mb] = GridGenerator(0,1,1,2,5,10,0);

        % ----------------------------------------------------------------
        %                      Concatenate grids
        % ----------------------------------------------------------------

        % Concatenate arrays of x coordinates of nodes
        xx = [xxa; xxb];

        % Concatenate arrays of y coordinates of nodes
        yy = [yya; yyb];

        % Concatenate arrays of x coordinates of cell vertices 
        xx1 = [xx1a; xx1b];

        % Concatenate arrays of y coordinates of cell vertices 
        yy1 = [yy1a; yy1b];

        % Concatenate arrays of cell areas 
        VV = [VVa; VVb];

        % Number of nodes
        Nnodes = length(xx);

        % Check consistency of grid spacings
        if dxa ~= dxb
            error('Grid spacings in x-direction are different.')
        else
            dx = dxa;
        end

        if dya ~= dyb
            error('Grid spacings in y-direction are different.')
        end

        % Find maximum vertex number in grid a
        maxMa = max(Ma(:));

        % Concatenate arrays that map each node to its cell vertices
        % (shift values of Mb array by maximum vertex number in grid a)
        M = [Ma; Mb + maxMa];

        % ----------------------------------------------------------------
        %                         Plot grid
        % ----------------------------------------------------------------

        % Change directory
        cd(Directory)

        PlotGrid(xx,yy,M,xx1,yy1,'b','k')

        % ----------------------------------------------------------------
        %                     Generate neighbor list
        % ----------------------------------------------------------------
       
        % Change directory
        cd ../Source/

        % Horizon
        del = 2.0*dx;

        [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList([],[],xx,yy,xx1,yy1,M,del,[],[],VV,omega,AlgName,flag_RDUG);

        % ----------------------------------------------------------------
        %              Break bonds crossing domain boundary
        % ----------------------------------------------------------------

        % Coordinates of endpoints of boundary edge representing a "prenoth" 
        Xc1 = 1;
        Yc1 = 1;
        Xc2 = 1;
        Yc2 = 2;

        u_NA = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);

    % --------------------------------------------------------------------
    %                         Circular domain
    % --------------------------------------------------------------------

    elseif strcmp(shape,'Circle')

        % ----------------------------------------------------------------
        %                 Generate square domain grid
        % ----------------------------------------------------------------

        % Change directory
        cd ../Source/

        % Generate uniform grid over a square domain
        [xx,yy,~,~,dx,~,VV,xx1,yy1,M] = GridGenerator(0,1,0,1,15,15,0);
        
        % Number of nodes
        Nnodes = length(xx);

        % ----------------------------------------------------------------
        %                  Remove nodes outside circle
        % ----------------------------------------------------------------

        % Circular domain center coordinates (Xc,Yc) and radius R
        Xc = 0.5;
        Yc = 0.5;
        R  = 0.5;

        % Initialize temporary arrays
        xxtemp = [];
        yytemp = [];
        VVtemp = [];
        Mtemp  = [];

        % Loop through all nodes ui
        for ui  = 1:Nnodes

            % Get x-coordinate of node ui
            xi = xx(ui);

            % Get y-coordinate of node ui
            yi = yy(ui);

            % Get area of cell corresponding to node ui
            VVi = VV(ui); 

            % Check if node ui is inside circular domain
            if (xi-Xc)^2 + (yi-Yc)^2 < R^2

                % Update temporary arrays
                xxtemp = [xxtemp; xi];     % x-coordinates of nodes inside circular domain
                yytemp = [yytemp; yi];     % y-coordinates of nodes inside circular domain
                VVtemp = [VVtemp; VVi];    % areas of cells corresponding to nodes inside circular domain
                Mtemp  = [Mtemp; M(ui,:)]; % array that maps each node inside circular domain to its cell vertices 
             
            end
        end

        % Set arrays with updated values
        xx = xxtemp;
        yy = yytemp;
        VV = VVtemp;
        M  = Mtemp;

        % Number of nodes
        Nnodes = length(xx);

        % ----------------------------------------------------------------
        %           Remove vertices that do not belong to cells 
        %              corresponding to nodes inside circle
        % ----------------------------------------------------------------

        % Find array of vertices' numbers corresponding to nodes inside circle
        I_vertices = unique(M(:));

        % Set cell vertices arrays with updated values
        xx1 = xx1(I_vertices);
        yy1 = yy1(I_vertices);

        % Loop over vertices' numbers
        for n = 1:length(I_vertices)

            % Vertex number
            nv = I_vertices(n);

            % Update array M: replace values of entries with "old" vertex number (nv) 
            %                 by "new" vertex number (n)
            I = (M==nv);
            M = M + (n - nv)*I;

        end

        % ----------------------------------------------------------------
        %                         Plot grid
        % ----------------------------------------------------------------

        % Change directory
        cd(Directory)

        % Plot grid
        PlotGrid(xx,yy,M,xx1,yy1,'b','k')

        % Box
        box('on')

        % ----------------------------------------------------------------
        %                     Generate neighbor list
        % ----------------------------------------------------------------

        % Change directory
        cd ../Source/

        % Horizon
        del = 2.0*dx;

        [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList([],[],xx,yy,xx1,yy1,M,del,[],[],VV,omega,AlgName,flag_RDUG);

    % --------------------------------------------------------------------
    %              Square with a circular hole at its center
    % --------------------------------------------------------------------

    elseif strcmp(shape,'SquareWithHole')

        % ----------------------------------------------------------------
        %                  Generate square domain grid
        % ----------------------------------------------------------------

        % Change directory
        cd ../Source/

        [xx,yy,~,~,dx,~,VV,xx1,yy1,M] = GridGenerator(0,1,0,1,20,20,0);

        % Number of nodes
        Nnodes = length(xx);

        % ----------------------------------------------------------------
        %                  Remove nodes inside hole
        % ----------------------------------------------------------------

        % Hole center coordinates (Xc,Yc) and radius R
        Xc = 0.5;
        Yc = 0.5;
        R  = 0.13;

        % Initialize temporary arrays
        xxtemp = [];
        yytemp = [];
        VVtemp = [];
        Mtemp  = [];

        % Loop through all nodes ui
        for ui  = 1:Nnodes

            % Get x-coordinate of node ui
            xi = xx(ui);

            % Get y-coordinate of node ui
            yi = yy(ui);

            % Get area of cell corresponding to node ui
            VVi = VV(ui); 

            % Check if node ui is outside (circular) hole
            if (xi-Xc)^2 + (yi-Yc)^2 > R^2

                % Update temporary arrays
                xxtemp = [xxtemp; xi];     % x-coordinates of nodes outside hole
                yytemp = [yytemp; yi];     % y-coordinates of nodes outside hole
                VVtemp = [VVtemp; VVi];    % areas of cells corresponding to nodes outside hole
                Mtemp  = [Mtemp; M(ui,:)]; % array that maps each node outside hole to its cell vertices
             
            end
        end

        % Set arrays with updated values
        xx = xxtemp;
        yy = yytemp;
        VV = VVtemp;
        M  = Mtemp;

        % Number of nodes
        Nnodes = length(xx);

        % ----------------------------------------------------------------
        %           Remove vertices that do not belong to cells 
        %              corresponding to nodes outside hole
        % ----------------------------------------------------------------

        % Find array of vertices' numbers corresponding to nodes outside hole
        I_vertices = unique(M(:));

        % Set cell vertices arrays with updated values
        xx1 = xx1(I_vertices);
        yy1 = yy1(I_vertices);

        % Loop over vertices' numbers
        for n = 1:length(I_vertices)

            % Vertex number
            nv = I_vertices(n);

            % Update array M: replace values of entries with "old" vertex number (nv) 
            %                 by "new" vertex number (n)
            I = (M==nv);
            M = M + (n - nv)*I;

        end

        % ----------------------------------------------------------------
        %                         Plot grid
        % ----------------------------------------------------------------

        % Change directory
        cd(Directory)

        % Plot grid
        PlotGrid(xx,yy,M,xx1,yy1,'b','k')

        % ----------------------------------------------------------------
        %                     Generate neighbor list
        % ----------------------------------------------------------------

        % Change directory
        cd ../Source/

        % Horizon
        del = 3.5*dx;

        [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList([],[],xx,yy,xx1,yy1,M,del,[],[],VV,omega,AlgName,flag_RDUG);

        % ----------------------------------------------------------------
        %                Break bonds intersecting hole
        % ----------------------------------------------------------------        
        
        % Note: we approximate the circular hole as a twenty-sided polygon
        %       for the purpose of breaking bonds

        % Number of polygonal sides
        Nsides = 20; 

        % Angle spacing for triangular decomposition of polygon
        dtheta = 2*pi/Nsides;

        % Loop through polygonal sides
        for n = 1:Nsides

            % Angle corresponding to one endpoint of polygonal side
            theta = (n-1)*dtheta;

            % Coordinates of one endpoint of polygonal side
            Xc1 = Xc + R*cos(theta);
            Yc1 = Yc + R*sin(theta);

            % Coordinates of the other endpoint of polygonal side
            Xc2 = Xc + R*cos(theta + dtheta);
            Yc2 = Yc + R*sin(theta + dtheta);

            % Break bonds intersecting polygonal side 
            u_NA = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);

            % Draw line representing polygonal side 
            line([Xc1 Xc2],[Yc1 Yc2],'LineStyle',':','LineWidth',2,'Color','r')

        end

    else

        error('Shape unknown.')

    end

    % --------------------------------------------------------------------
    %                          Plot bonds
    % --------------------------------------------------------------------

    if flag_PlotBonds == 1

        % Find maximum number of neighbors any node can have
        zmax = length(u_NA(1,:));

        % Loop through all nodes ui
        for ui = 1:Nnodes

            % Get x-coordinate of node ui
            xi = xx(ui);

            % Get y-coordinate of node ui
            yi = yy(ui);

            % Loop over neighbors of node ui
            for z = 1:zmax

                % Get neighbor cell number
                uk = u_NA(ui,z);

                % Only consider neighbor cells
                if uk ~= 0

                    % Get x-coordinate of quadrature point
                    xk = x_hat_NA(ui,z);

                    % Get y-coordinate of quadrature point
                    yk = y_hat_NA(ui,z);

                    % Draw line representing a bond
                    line([xi xk],[yi yk],'Color','cyan','LineStyle','-')

                end

            end
        end

    end

    % Change directory
    cd(Directory)
    
end
