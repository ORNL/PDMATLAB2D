
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
% The function PlotNeighborAreasCoords plots the family of a node at the 
% origin (source node i) and the associated cell areas
% ========================================================================

% Input
% -----
% m       : m-ratio (horizon divided by grid spacing)
% AlgName : algorithm name in string format
%           'FA'    : FA algorithm
%           'PA-AC' : PA-AC algorithm
%           'IPA-AC': IPA-AC algorithm

function PlotNeighborAreasCoords(m,AlgName)
    
    % Identify current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % Tolerance
    tol = 1E-15;

    % Source node i coordinates (at the origin)
    xi = 0;
    yi = 0;

    % Horizon
    del = 1.0;

    % ----------------------------------------------------------------
    %                        Generate grid
    % ----------------------------------------------------------------

    % Number of nodes per dimension
    Nx = 2*ceil(m) + 1;
    Ny = Nx;

    % Uniform grid spacing
    h = del/m;

    % Domain limits in the x-direction
    Xo = xi - ceil(m)*h - 0.5*h;
    Xn = xi + ceil(m)*h + 0.5*h;

    % Domain limits in the y-direction
    Yo = yi - ceil(m)*h - 0.5*h;
    Yn = yi + ceil(m)*h + 0.5*h;

    % Grid perturbation coefficient
    PG = 0;

    % Generate grid
    [xx,yy,~,~,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);

    % ----------------------------------------------------------------
    %                          Plot grid 
    % ----------------------------------------------------------------

    % Change directory
    cd(Directory)

    PlotGrid(xx,yy,M,xx1,yy1,'k','k')
    hold on
    axis([Xo Xn Yo Yn])

    % ----------------------------------------------------------------
    %  Plot neighbor areas and centroids (in case of IPA-AC algorithm)
    % ----------------------------------------------------------------

    % Change directory
    cd ../Source/

    % Number of cells
    Ncells = length(xx);

    % Cell area color (gray)
    carea = [0.75 0.75 0.75];

    % Run over cells
    for k = 1:Ncells

        % Coordinates of node k
        xk = xx(k);
        yk = yy(k);

        % Compute distance squared between source node i and node k
        r2 = ((xk-xi)^2+(yk-yi)^2);

        % Check if node k is different from source node i
        if r2 > tol

            % Coordinates of four vertices of cell k
            xx1k = xx1(M(k,:));
            yy1k = yy1(M(k,:));

            % Volume of cell k
            VVk = VV(k);

            % Compute neighbor area and coordinates of quadrature point associated with the cell k relative to source node i
            [Vk,~,xk_hat,yk_hat] = NeighborAreaBondLengthCoord(xi,xk,yi,yk,dx,dy,xx1k,yy1k,del,VVk,AlgName);

            % Check if cell k is associated with the neighborhood of source node i
            if Vk > 0

                % Coordinates of four vertices of cell k
                xvertex = xx1k;
                yvertex = yy1k;

                if strcmp(AlgName,'FA')  

                    % Plot cell k area
                    patch(xvertex,yvertex,'','FaceColor',carea)

                    % Plot node k again for better visualization
                    plot(xk,yk,'ob','MarkerFaceColor','b')

                else

                    % Check if cell k is inside neighborhood of source node i
                    if abs(Vk - h^2) < tol

                        % Plot cell k area
                        patch(xvertex,yvertex,'','FaceColor',carea)

                        % Plot node k again for better visualization
                        plot(xk,yk,'ob','MarkerFaceColor','b')

                        if strcmp(AlgName,'IPA-AC')     
                            % Plot neighbor centroid
                            plot(xk_hat,yk_hat,'ob','MarkerFaceColor','w')
                        end

                    else

                        % ------------------------------------
                        % Use Quasi Monte Carlo (QMC) Approach
                        % ------------------------------------

                        % Number of QMC points
                        N_QMC = 1E+5;

                        % QMC random points
                        Q = qrandstream('sobol',2);
                        p = qrand(Q,N_QMC);

                        % Map points inside cell k
                        minxk = min(xvertex);
                        minyk = min(yvertex);

                        x_QMC = minxk + dx.*p(:,1);
                        y_QMC = minyk + dy.*p(:,2);

                        % Distance squared between source node i and QMC points
                        r2_QMC = (x_QMC - xi).^2 + (y_QMC - yi).^2;

                        % Find QMC points inside neighborhood of source node i
                        mask = find( r2_QMC < del^2 );

                        % Plot QMC points inside neighborhood of source node i
                        plot(x_QMC(mask),y_QMC(mask),'.','color',carea)

                        % Plot node k again for better visualization
                        plot(xk,yk,'ob','MarkerFaceColor','b')

                        if strcmp(AlgName,'IPA-AC')     
                            % Plot neighbor centroid
                            plot(xk_hat,yk_hat,'ob','MarkerFaceColor','w')
                        end
                    end

                end

            end

        end

    end

    % Plot cell k edges again for better visualization
    for i = 1:3
        line([xx1(M(:,i)) xx1(M(:,i+1))]', [yy1(M(:,i)) yy1(M(:,i+1))]','Color','k')
    end
    line([xx1(M(:,4)) xx1(M(:,1))]', [yy1(M(:,4)) yy1(M(:,1))]','Color','k')

    % Plot boundary of neighborhood of source node i
    theta = 0:0.001:2*pi;
    xboundary = del*cos(theta);
    yboundary = del*sin(theta);

    plot(xboundary,yboundary,'-b')

    % Change directory
    cd(Directory)

end
