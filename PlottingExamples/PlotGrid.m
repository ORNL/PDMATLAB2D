
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
% The function PlotGrid plots the grid
% ========================================================================

% Input
% -----
% xx     : x coordinates of all nodes in the grid (1D array of length Nx x Ny)
% yy     : y coordinates of all nodes in the grid (1D array of length Nx x Ny)
% M      : array that maps each node to its cell vertices 
%          (2D array of size (Nx x Ny)-by-4)
% xx1    : x coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% yy1    : y coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% cnodes : color for nodes 
% cedges : color for cell edges 

function PlotGrid(xx,yy,M,xx1,yy1,cnodes,cedges)

    % Create figure
    figure1 = figure;
    axes1 = axes('Parent',figure1);

    % Plot nodes
    scatter(xx,yy,cnodes,'filled')

    % Plot cell edges
    for i = 1:3
        line([xx1(M(:,i)) xx1(M(:,i+1))]', [yy1(M(:,i)) yy1(M(:,i+1))]','Color',cedges)
    end
    line([xx1(M(:,4)) xx1(M(:,1))]', [yy1(M(:,4)) yy1(M(:,1))]','Color',cedges)

    % Set axes font size and latex font style
    set(axes1,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex')

    % Axes labels
    xlabel('$x$','Interpreter','latex','FontSize',30)
    ylabel('$y$','Interpreter','latex','FontSize',30)

    % Find domain limits
    minx = min(xx1);
    maxx = max(xx1);
    miny = min(yy1);
    maxy = max(yy1);

    % Use real aspect ratio
    pbaspect([(maxx - minx) (maxy - miny) 1])

end
