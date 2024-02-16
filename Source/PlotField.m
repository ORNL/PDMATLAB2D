
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
% The function PlotField plots the nodes colored by a given field
% ========================================================================

% Input
% -----
% nfig    : figure number
% xx      : x coordinates of all nodes in the grid 
% yy      : y coordinates of all nodes in the grid 
% cnodes  : field values of all nodes in the grid (for color assignment)
% ctitle  : colorbar title
% psize   : point size
% climits : 1D array containing the color limits as [cmin cmax]
% cmap    : colormap
% box     : 1D array containing the box limits as [xmin xmax ymin ymax]

% Discussion
% ----------
% The function can take empty arrays for the following inputs:
% box     : if empty, then the axes limits are set to [min(xx) max(xx) min(yy) max(yy)]
% cmap    : if empty, then the function uses MATLAB's default colormap
% climits : if empty, then the color limits are set to [min(cnodes)-tol max(cnodes)+tol] 
%           where tol is a tolerance

function PlotField(nfig,xx,yy,cnodes,ctitle,psize,climits,cmap,box)

    % Tolerance
    tol = 1E-15;

    % Create figure
    fig   = figure(nfig);
    clf(fig)                     % clear current figure
    set(fig,'color','w');        % use white background
    axfig = axes('Parent',fig);

    % Box limits
    if isempty(box)
        minx = min(xx); maxx = max(xx);
        miny = min(yy); maxy = max(yy);
    else
        minx = box(1); maxx = box(2);
        miny = box(3); maxy = box(4);
    end

    % Plot nodes colored by field values
    scatter(xx,yy,psize,cnodes,'filled')

    % Set axes font size and latex font style
    set(axfig,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex')

    % Axes labels
    xlabel('$x$','Interpreter','latex','FontSize',30)
    ylabel('$y$','Interpreter','latex','FontSize',30)

    % Use real aspect ratio
    pbaspect([(maxx - minx) (maxy - miny) 1])

    % Axes limits
    axis([minx maxx miny maxy])

    % Colorbar
    hcb = colorbar;

    % Colorbar title latex style
    hcb.Title.Interpreter = 'latex';

    % Colorbar ticks latex style
    set(hcb,'TickLabelInterpreter','latex')

    % Colorbar title
    hcb.Title.String = ctitle;

    % Colorbar title font size
    hcb.Title.FontSize = 20;

    % Get MATLAB version year
    MATLABversion = version('-release'); 
    MATLAByear    = str2num(MATLABversion(1:4));
    
    % Color limits
    
    % Note: the MATLAB function 'caxis' was renamed 'clim' starting with
    %       version R2022a (https://www.mathworks.com/help/matlab/ref/clim.html)

    if MATLAByear < 2022
        % Use caxis to set colormap limits
        if isempty(climits)
            caxis([min(cnodes)-tol max(cnodes)+tol])
        else
            caxis(climits)
        end
    else
        % Use clim to set colormap limits
        if isempty(climits)
            clim([min(cnodes)-tol max(cnodes)+tol])
        else
            clim(climits)
        end
    end

    % Colormap
    if isempty(cmap)
        colormap;
    else
        colormap(cmap);
    end 

    % Ticks label format
    hcb.Ruler.TickLabelFormat = '%.1f';

end

