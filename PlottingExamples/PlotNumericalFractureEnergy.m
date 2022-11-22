
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
% The function PlotNumericalFractureEnergy plots the setting for the 
% numerical computation of the fracture energy
% ========================================================================

% Input
% -----
% del : horizon
% m   : m-ratio (horizon divided by grid spacing)

function PlotNumericalFractureEnergy(del,m)

    % Current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                     Domain and discretization
    % --------------------------------------------------------------------
    
    % Uniform grid spacing
    h = del/m;

    % Coordinates of intersection point between horizontal and vertical lines
    xI = 0;
    yI = 0;

    % Domain limits
    minx = -1;
    maxx =  1;
    miny = -1-0.5*h;
    maxy =  1+0.5*h;

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
    % Create "intersection boundary": intersection between neighborhood 
    %  boundary of rightmost line node and right side of vertical line
    % -------------------------------------------------------------------- 

    % Full neighborhood boundary for node (xI-0.5*dx, yI)
    theta = 0:pi/1000:2*pi;
    xbdy  = xI - 0.5*dx + del*cos(theta);
    ybdy  = yI + del*sin(theta);
    
    % Intersection between neighborhood boundary and right side
    ybdy  = ybdy(xbdy > xI - 1e-14);
    xbdy  = xbdy(xbdy > xI - 1e-14);

    % --------------------------------------------------------------------
    %                         Plot the system 
    % -------------------------------------------------------------------- 
    
    % Create figure
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');

    % Plot vertical line (vertical surface)
    line([xI xI],[miny maxy],'Color',[255 0 255]/255,'LineStyle','-','LineWidth',4)

    % Plot left-side line with nodes
    line([xI-del xI],[yI yI],'Color','k')
    plot(xxL,yyL,'og','MarkerFaceColor','g') 

    % Plot "intersection boundary" and corresponding radii
    plot(xbdy,ybdy,'.b','LineWidth',1)
    H = sqrt(del^2 - (0.5*dx)^2);
    line([xI-0.5*dx xI],[yI yI-H],'Color','k','LineStyle','-')
    line([xI-0.5*dx xI],[yI yI+H],'Color','k','LineStyle','-')

    % Plot right-side grid nodes and corresponding cells
    plot(xxR,yyR,'ok','MarkerFaceColor','k')
    
    for j = 1:dimR

        % Node coordinates
        x = xxR(j);
        y = yyR(j);

        % Left cell edge
        line([x-0.5*h x-0.5*h],[y-0.5*h y+0.5*h],'Color',[17 17 17]/255,'LineStyle','-','LineWidth',1)

        % Right cell edge
        line([x+0.5*h x+0.5*h],[y-0.5*h y+0.5*h],'Color',[17 17 17]/255,'LineStyle','-','LineWidth',1)

        % Top cell edge
        line([x-0.5*h x+0.5*h],[y+0.5*h y+0.5*h],'Color',[17 17 17]/255,'LineStyle','-','LineWidth',1)

        % Bottom cell edge
        line([x-0.5*h x+0.5*h],[y-0.5*h y-0.5*h],'Color',[17 17 17]/255,'LineStyle','-','LineWidth',1)
    end

    % Set axes font size and latex font style
    set(axes1,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex')

    % Axes labels
    xlabel('$x$','Interpreter','latex','FontSize',30)
    ylabel('$y$','Interpreter','latex','FontSize',30)

    % Box
    box('on')

    % Axes limits
    axis([minx maxx miny maxy])

    % Aspect ratio
    pbaspect([2 2+h 1])

    % Change directory
    cd(Directory)

end