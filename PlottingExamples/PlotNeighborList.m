
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
% The function PlotNeighborList plots the bonds for all nodes in a domain 
% and emphasizes a single neighborhood
% ========================================================================

% Input
% -----
% ncase : case number  
%         1: Rectangular domain, uniform grid with FA algorithm
%         2: Rectangular domain, uniform grid with PA-AC algorithm
%         3: Rectangular domain, uniform grid with IPA-AC algorithm
%         4: Rectangular domain, regular (nonuniform) grid with FA algorithm
%         5: Rectangular domain, perturbed regular grid with FA algorithm
%         6: Circular domain, irregular grid with FA algorithm

function PlotNeighborList(ncase)

    % Identify current directory
    Directory = pwd;

    % Change directory
    cd ../Source/
    
    % Cases 1-3: Uniform grid
    if ncase == 1 || ncase == 2 || ncase == 3

        % Number of nodes per dimension
        Nx = 14;
        Ny = 7;

        % Grid perturbation coefficient
        PG = 0;

        % Flag for Rectangular Domain Uniform Grid (RDUG)
        flag_RDUG = 1;

    % Cases 4-5: Nonuniform regular grid (unperturbed (ncase = 4) or perturbed (ncase = 5))  
    elseif ncase == 4 || ncase == 5

        % Number of nodes per dimension
        Nx = 10;
        Ny = 10;
        
        % Grid perturbation coefficient
        if ncase == 4 
            PG = 0;
        else
            PG = 0.99;
        end
        
        % Flag for Rectangular Domain Uniform Grid (RDUG)
        flag_RDUG = 0;

    elseif ncase == 6

        % Flag for Rectangular Domain Uniform Grid (RDUG)
        flag_RDUG = 0;

    else

        % Change directory
        cd(Directory)

        error('Unknown case number.')

    end

    % Algorithm for computation of neighbor areas
    if ncase == 2 
        AlgName = 'PA-AC';
    elseif ncase == 3 
        AlgName = 'IPA-AC';
    else
        AlgName = 'FA';
    end

    if ncase == 6

        % ----------------------------------------------------------------
        %        Generate irregular grid over circular domain
        % ----------------------------------------------------------------

        % Choose circular domain center coordinates (Xc,Yc) and radius R
        Xc = 1;
        Yc = 0.5;
        R  = 1;

        % Estimated number of vertices
        Nvertices = 70;

        % Generate grid
        [xx,yy,VV,xx1,yy1,M] = CircularDomainQRandomGrid(Xc,Yc,R,Nvertices);

        % Horizon
        del = 0.4;

        % Domain center coordinates
        midx = Xc;
        midy = Yc;

        % Number of nodes
        Nnodes = length(xx);

        % ----------------------------------------------------------------
        %                          Plot grid
        % ----------------------------------------------------------------

        % Create figure
        figure1 = figure;
        axes1 = axes('Parent',figure1);
        hold all

        % Plot nodes
        scatter(xx,yy,'b','filled')

        % Define domain boundary
        theta = 0:0.01:2*pi;
        xbdy = Xc + R*cos(theta);
        ybdy = Yc + R*sin(theta);

        % Plot domain boundary
        plot(xbdy,ybdy,'-k','LineWidth',1.15)

        % Plot cells
        triplot(M,xx1,yy1,'k');

        % Set axes font size and latex font style
        set(axes1,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex')

        % Axes labels
        xlabel('$x$','Interpreter','latex','FontSize',30)
        ylabel('$y$','Interpreter','latex','FontSize',30)

        % Set axes limits
        axis([midx-R midx+R midy-R midy+R])

        % Use real aspect ratio
        pbaspect([1 1 1])

        % Add box
        box('on')

        % ----------------------------------------------------------------
        %                     Generate neighbor list
        % ----------------------------------------------------------------

        % Influence function order indicator
        omega = 0;

        [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList([],[],xx,yy,xx1,yy1,M,del,[],[],VV,omega,AlgName,flag_RDUG);

    else
        
        % ----------------------------------------------------------------
        %           Generate (possibly perturbed) regular grid 
        %                over a rectangular domain
        % ----------------------------------------------------------------

        % Domain limits in the x-direction
        Xo = 0;
        Xn = 2;

        % Domain limits in the y-direction
        Yo = 0;
        Yn = 1;

        % Generate grid
        [xx,yy,~,~,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);

        % Number of nodes
        Nnodes = Nx*Ny;

        % Cases 1-3: Uniform grid
        if ncase == 1 || ncase == 2 || ncase == 3
            del = 2.5*dx;

        % Cases 4-5: Nonuniform regular grid (unperturbed (ncase = 4) or perturbed (ncase = 5))
        else
            del = 3.5*dy;

        end

        % Domain center coordinates
        midx = (Xn - Xo)/2;
        midy = (Yn - Yo)/2;

        % ----------------------------------------------------------------
        %                         Plot grid
        % ----------------------------------------------------------------

        % Change directory
        cd(Directory)

        PlotGrid(xx,yy,M,xx1,yy1,'b','k')
        hold on
        axis([Xo Xn Yo Yn])

        % ----------------------------------------------------------------
        %                    Generate neighbor list
        % ----------------------------------------------------------------

        % Change directory
        cd ../Source/

        % Influence function order indicator
        omega = 0;

        [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);

    end

    % ----------------------------------------------------------------
    %            Find node nearest to the domain center
    % ----------------------------------------------------------------

    % Find distance squared of all nodes to the domain center
    D2 = abs(xx - midx).^2 + abs(yy - midy).^2;

    % Choose first node among closest to the domain center
    [~, Isort] = sort(D2);
    umid = Isort(1);

    % ----------------------------------------------------------------
    %                       Plot all bonds 
    % ----------------------------------------------------------------

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

            % Get neighbor node number
            uk = u_NA(ui,z);
            
            if uk ~= 0

                % Get x-coordinate of quadrature point
                xk = x_hat_NA(ui,z);
                
                % Get y-coordinate of quadrature point
                yk = y_hat_NA(ui,z);

                % Draw line representing a bond
                line([xi xk],[yi yk],'Color','cyan','LineStyle','-')

                hold on

            end

        end
    end

    % Plot nodes again for better visualization
    scatter(xx,yy,'b','filled')

    % ----------------------------------------------------------------
    %           Highlight bonds and neighborhood boundary 
    %       of a source node among the closest to the domain center
    % ----------------------------------------------------------------

    % Set node number
    ui = umid;

    % Get x-coordinate of node ui
    xi = xx(ui);
    
    % Get y-coordinate of node ui
    yi = yy(ui);

    % Loop over neighbors of node ui
    for z = 1:zmax

        % Get neighbor node number
        uk = u_NA(ui,z);

        if uk ~= 0

            % Get x-coordinate of quadrature point
            xk = x_hat_NA(ui,z);

            % Get y-coordinate of quadrature point
            yk = y_hat_NA(ui,z);

            % Draw line representing a bond
            line([xi xk],[yi yk],'Color','red','LineStyle','-','LineWidth',1.15)

            % Plot quadrature point
            scatter(xk,yk,'r','filled')

            hold on

        end

    end

    % Plot source node again for better visualization
    scatter(xi,yi,'k','filled')

    % Define neighborhood boundary
    theta = 0:0.01:2*pi;
    xbdy = xi + del*cos(theta);
    ybdy = yi + del*sin(theta);

    % Plot neighborhood boundary
    plot(xbdy,ybdy,'-r','LineWidth',1.15)

    % Change directory
    cd(Directory)

end

% ========================================================================
% The function CircularDomainQRandomGrid generates a quasi-random grid of
% nodes with corresponding areas within a circular domain
% ========================================================================

% Input
% -----
% Xc        : x-coordinate of center of domain
% Yc        : y-coordinate of center of domain
% R         : domain radius
% Nvertices : estimation of number of vertices for triangulation

% Output 
% ------
% xx        : array of x coordinates of all nodes in the grid (1D array)
% yy        : array of y coordinates of all nodes in the grid (1D array)
% VV        : array of areas of cells corresponding to the nodes in the grid (1D array)
% xx1       : array of x coordinates of all cell vertices (1D array)
% yy1       : array of y coordinates of all cell vertices (1D array)
% M         : array that maps each node to its cell vertices 
%             (2D array of size (number of nodes)-by-3)

% Discussion:
% ----------
% Nvertices provides only an estimation of the number of cell vertices 
% in the grid. The actual number of cell vertices is generated with the 
% following steps:
%
% Step 1: We generate quasi-random points using a Quasi Monte Carlo (QMC) 
%         approach within a circumscribed square of the circular domain, 
%         i.e., a square of edge length 2R centered at the center of the domain. 
%         The number of QMC points in the square is approximated as follows:
%
%         Nvertices * (Area of square)/(Area of circle)
%
% Step 2: We select all QMC points within a distance of 0.9R to the
%         center of the domain as cell vertices. The choice of 0.9R is to 
%         select vertices not too close to the domain boundary (we select 
%         vertices on the domain boundary in Step 3).
%
% Step 3: We select vertices on the domain boundary. The number of such 
%         vertices is computed assuming the domain boundary is "represented" 
%         by an annulus of width 0.1R adjacent to the domain boundary (and 
%         inside the domain) and using the following approximation:
%
%         Nvertices * (Area of annulus)/(Area of circle) ~ 0.19 * Nvertices

function [xx,yy,VV,xx1,yy1,M] = CircularDomainQRandomGrid(Xc,Yc,R,Nvertices)

        % ----------------------------------------------------------------
        %                   Generate cell vertices
        % ----------------------------------------------------------------

        % Number of QMC points in circumscribed square 
        N_QMC = floor(Nvertices*((2*R)^2)/(pi*R^2));

        % QMC random points
        Q = qrandstream('sobol',2);
        p = qrand(Q,N_QMC);

        % Map QMC points inside circumscribed square
        minxc = Xc - R;
        minyc = Yc - R;

        x_QMC = minxc + 2*R.*p(:,1);
        y_QMC = minyc + 2*R.*p(:,2);

        % Distance squared between domain center and QMC points
        r2_QMC = (x_QMC - Xc).^2 + (y_QMC - Yc).^2;

        % Find QMC points inside "inner domain" (circle of radius 0.9R centered
        % at the domain center)
        mask = find( r2_QMC < (0.9*R)^2 );

        xx1 = x_QMC(mask);
        yy1 = y_QMC(mask);

        % Estimate number of vertices on domain boundary
        NverticesBdy = floor(0.19*Nvertices);

        % Generate vertices on domain boundary
        dtheta = 2*pi/NverticesBdy;
        theta = 0:dtheta:2*pi;
        xbdy = Xc + R*cos(theta);
        ybdy = Yc + R*sin(theta);

        % Add vertices on domain boundary to array of vertices inside "inner domain"
        xx1 = [xx1' xbdy]';
        yy1 = [yy1' ybdy]';

        % ----------------------------------------------------------------
        %                Perform Delaunay triangulation
        % ----------------------------------------------------------------

        M = delaunay(xx1,yy1);

        % ----------------------------------------------------------------
        %       Compute areas and centroids of triangular cells
        % ----------------------------------------------------------------

        % Size of M
        [s1,~] = size(M);
        Nnodes = s1; % Note: number of cells = number of nodes

        % Initialze arrays of x- and y- coordinates as well as array of cell areas of all nodes in the grid
        xx = zeros(Nnodes,1);
        yy = zeros(Nnodes,1);
        VV = zeros(Nnodes,1);

        % Loop over cells
        for ui = 1:Nnodes

            % Coordinates of three vertices of cell ui
            xx1i = xx1(M(ui,:));
            yy1i = yy1(M(ui,:));

            % Centroid coordinates
            xx(ui) = sum(xx1i)/3;
            yy(ui) = sum(yy1i)/3;

            % Cell area
            VV(ui) = polyarea(xx1i,yy1i);

        end

end
