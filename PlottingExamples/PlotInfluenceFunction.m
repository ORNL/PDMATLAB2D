
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
% The function PlotInfluenceFunction plots the influence functions
% ========================================================================

% Input
% -----
% del               : horizon
% rmax              : maximum r
% dr_fine           : plot resolution for lines
% dr_coarse         : plot resolution for markers
% omega_array_input : vector of influence function order indicators to plot

function PlotInfluenceFunction(del,rmax,dr_fine,dr_coarse,omega_array_input)

    % Plotting functions directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % Tolerance
    tol = 1E-15;

    % Note: We divide the array of r values into r <= delta and r > delta
    %       for more accurate visualization of the piecewise constant 
    %       influence function.

    % Define coarse vector for markers
    r_array_markers_l = 0:dr_coarse:del;                                   % Points with r <= delta
    r_array_markers_r = r_array_markers_l(end) + dr_coarse:dr_coarse:rmax; % Points with r > delta

    % Define fine vectors for lines
    r_array_lines_l = 0:dr_fine:del;        % Points with r <= delta
    r_array_lines_r = del+tol:dr_fine:rmax; % Points with r > delta

    % Create figure
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold all

    % Vector of influence function order indicators
    omega_array = [0 0.5 1 3 5 7];
    
    % Vector of markers
    marker_array = ['x' 'o' '<' '>' '+' '*'];

    % Vector of colors
    color_array = [  0    0    0; % black
                     0    0    1; % blue
                     1    0    0; % red
                     1    0    1; % magenta
                     0 128/255 0; % green
                     1 140/255 0; % dark orange
                   ];

    % Vector of legends
    legend_array = {'$\omega_0$','$\omega_{0.5}$','$\omega_{1}$','$\omega_{3}$','$\omega_{5}$','$\omega_{7}$'};

    % Initialize legend array for inputs
    legend_array_input = [];

    % --------------------------------------------------------------------
    %                          Plot markers
    % --------------------------------------------------------------------

    % ------------------------
    %       r <= delta
    % ------------------------

    % Note: We plot markers using a coarse grid, instead of using the fine
    %       grid, so that the plot is more clear visually.
    %
    %       To obtain markers with lines in the plot legend, we then
    %       proceed as follows:
    %       (1) we first plot the markers with lines using the coarse grid,
    %       (2) we then "delete" the lines by plotting white lines on top.

    for i = 1:length(omega_array_input)

        % Find index of influence function to plot
        n = find(omega_array == omega_array_input(i));

        if ~isempty(n)

            % Evaluate influence function
            omega = InfluenceFunction(omega_array(n),r_array_markers_l,del);

            % Plot influence function (lines + markers)
            plot(r_array_markers_l,omega,marker_array(n),'LineWidth',0.1,'LineStyle','-','Color',color_array(n,:),'MarkerSize',8);

            % Assign legend
            legend_array_input = [legend_array_input legend_array(n)];

        end

    end

    for i = 1:length(omega_array_input)

        % Find index of influence function to plot
        n = find(omega_array == omega_array_input(i));

        if ~isempty(n)

            % Evaluate influence function
            omega = InfluenceFunction(omega_array(n),r_array_markers_l,del);

            % Plot influence function (white lines) 
            plot(r_array_markers_l,omega,'LineWidth',1,'LineStyle','-','Color','w');
        end
    end

    % ------------------------
    %        r > delta
    % ------------------------

    for i = 1:length(omega_array_input)

        % Find index of influence function to plot
        n = find(omega_array == omega_array_input(i));

        if ~isempty(n)

            % Evaluate influence function
            omega = InfluenceFunction(omega_array(n),r_array_markers_r,del);

            % Plot influence function (markers)
            plot(r_array_markers_r,omega,marker_array(n),'Color',color_array(n,:),'MarkerSize',8);

        end

    end

    % --------------------------------------------------------------------
    %                          Plot lines
    % --------------------------------------------------------------------

    for i = 1:length(omega_array_input)

        % Find index of influence function to plot
        n = find(omega_array == omega_array_input(i));

        if ~isempty(n)

            % Evaluate influence function: r <= delta
            omega_l = InfluenceFunction(omega_array(n),r_array_lines_l,del);

            % Evaluate influence function: r > delta
            omega_r = InfluenceFunction(omega_array(n),r_array_lines_r,del);

            % Plot influence function (lines)
            plot(r_array_lines_l,omega_l,r_array_lines_r,omega_r,'-','Color',color_array(n,:))

        end

    end

    % Legend
    hleg = legend([legend_array_input,'']);
    set(hleg,'Interpreter','latex','FontSize',20)

    % Set axes font size and latex font style
    set(axes1,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex')

    % Axes labels
    xlabel('$r$','Interpreter','latex','FontSize',30)
    ylabel('$\omega(r)$','Interpreter','latex','FontSize',30)

    % Add box
    box('on')

    % Change directory
    cd(Directory)

end
