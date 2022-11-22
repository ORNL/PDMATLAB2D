
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
% The function PlotPreNotch demonstrates the creation of pre-notches by 
% creating pre-notches and plotting intact bonds and broken bonds in 
% different colors
% ========================================================================

% Input
% -----
% ncase                 : case number  
%                         1: Horizontal pre-notch
%                         2: Inclined pre-notch
%                         3: Multiple random pre-notches
% flag_plot_broken_bonds: if == 1, then the function plots broken bonds
%                         in addition to intact bonds
%                         if == 0, then the function does not plot broken bonds

function PlotPreNotch(ncase,flag_plot_broken_bonds)

    % Identify current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                      Generate and plot grid 
    % --------------------------------------------------------------------

    % Domain limits in the x-direction
    Xo = 0;
    Xn = 2;

    % Domain limits in the y-direction
    Yo = 0;
    Yn = 1;

    % Number of nodes per dimension
    Nx = 20;
    Ny = 10;

    % Grid perturbation coefficient
    PG = 0;

    % Generate uniform grid
    [xx,yy,~,~,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);

    % Change directory
    cd(Directory)

    % Plot grid
    PlotGrid(xx,yy,M,xx1,yy1,'b','k')

    % --------------------------------------------------------------------
    %                      Generate neighbor list
    % --------------------------------------------------------------------

    % Change directory
    cd ../Source/

    % Horizon
    del = 2.5*dx;

    % Influence function order indicator
    omega = 0;

    % Algorithm for computation of neighbor areas
    AlgName = 'FA';

    % Flag for Rectangular Domain Uniform Grid (RDUG)
    flag_RDUG = 1;

    [u_NA,~,~,~,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);

    % --------------------------------------------------------------------
    %     Create pre-notch(es) and plot corresponding line segment(s)
    % --------------------------------------------------------------------

    switch ncase
        
        case 1
            % ------------------------------------------------------------
            %                    Horizontal pre-notch
            % ------------------------------------------------------------

            % Coordinates of one endpoint of the pre-notch
            Xc1 = Xo;
            Yc1 = 0.5*(Yo + Yn);

            % Coordinates of the other endpoint of the pre-notch 
            Xc2 = 0.5*(Xo + Xn);
            Yc2 = Yc1;

            % Create pre-notch
            [u_NA_pn] = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);

            % Draw pre-notch line segment
            line([Xc1 Xc2],[Yc1 Yc2],'LineWidth',3,'Color','k')

        case 2
            % ------------------------------------------------------------
            %                   Inclined pre-notch
            % ------------------------------------------------------------

            % Coordinates of one endpoint of the pre-notch
            Xc1 = 0.45;
            Yc1 = 0.28;

            % Coordinates of the other endpoint of the pre-notch 
            Xc2 = 1.45;
            Yc2 = 0.78;

            % Create pre-notch
            [u_NA_pn] = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);

            % Draw pre-notch line segment
            line([Xc1 Xc2],[Yc1 Yc2],'LineWidth',3,'Color','k')

        case 3
            % ------------------------------------------------------------
            %               Multiple random pre-notches
            % ------------------------------------------------------------

            % Number of pre-notches
            Npn = 15;

            % Initialize updated array of neighbor numbers
            u_NA_pn = u_NA;

            % Minimum possible pre-notch length
            minl = 0.1;

            % Maximim possible pre-notch length
            maxl = 0.25;

            % Set seed for random number generator
            rng('default');

            % Run over pre-notches
            for n = 1:Npn

                % Coordinates of one endpoint of the pre-notch (random)
                % Note: the addition and subtraction of "maxl" ensures
                %       the prenotch does not extend outside the domain
                Xc1 = (Xo+maxl) + ((Xn-maxl) - (Xo+maxl))*rand(1);
                Yc1 = (Yo+maxl) + ((Yn-maxl) - (Yo+maxl))*rand(1);

                % Random pre-notch length
                l = minl + (maxl-minl)*rand(1); 

                % Random orientation
                theta = (2*pi)*rand(1);

                % Coordinates of the other endpoint of the pre-notch
                Xc2 = Xc1 + l*cos(theta);
                Yc2 = Yc1 + l*sin(theta);

                % Create pre-notch
                [u_NA_pn] = PreNotch(xx,yy,u_NA_pn,Xc1,Yc1,Xc2,Yc2);

                % Draw pre-notch line segment
                line([Xc1 Xc2],[Yc1 Yc2],'LineWidth',2,'Color','k')

            end

        otherwise

            % Change directory
            cd(Directory)

            error('Invalid pre-notch case number.');

    end

    % --------------------------------------------------------------------
    %                        Plot intact bonds 
    % --------------------------------------------------------------------

    % Number of nodes
    Nnodes = length(xx);

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

            % Get updated neighbor cell number (after pre-notch)
            uk_pn = u_NA_pn(ui,z);

            if uk_pn ~= 0

                % Get x-coordinate of quadrature point
                xk = x_hat_NA(ui,z);

                % Get y-coordinate of quadrature point
                yk = y_hat_NA(ui,z);

                % Draw line representing a bond
                line([xi xk],[yi yk],'Color','cyan','LineStyle','-')

            end

        end
    end

    if flag_plot_broken_bonds == 1

        % ----------------------------------------------------------------
        %                      Plot broken bonds
        % ----------------------------------------------------------------

        % Loop through all nodes ui
        for ui = 1:Nnodes

            % Get x-coordinate of node ui
            xi = xx(ui);

            % Get y-coordinate of node ui
            yi = yy(ui);

            % Loop over neighbors of node ui
            for z = 1:zmax

                % Get neighbor cell number (before pre-notch)
                uk = u_NA(ui,z);

                if uk ~= 0

                    % Get x-coordinate of quadrature point
                    xk = x_hat_NA(ui,z);

                    % Get y-coordinate of quadrature point
                    yk = y_hat_NA(ui,z);

                    % Get updated neighbor cell number (after pre-notch)
                    uk_pn = u_NA_pn(ui,z);

                    if uk_pn == 0

                        % Draw line representing a bond
                        line([xi xk],[yi yk],'Color','red','LineStyle','-')

                    end

                end

            end
        end

    elseif flag_plot_broken_bonds == 0

    else

        % Change directory
        cd(Directory)

        error('flag_plot_broken_bonds should be 0 or 1.')

    end

    % Change directory
    cd(Directory)

end
