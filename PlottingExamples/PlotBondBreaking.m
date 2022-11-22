
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
% The function PlotBondBreaking plots a square domain in the undeformed or 
% deformed configuration having a single horizontal layer of nodes being 
% deformed with or without introducing a no-fail zone
% ========================================================================

% Input
% -----
% flag_deformation : if == 1, then a linear deformation in the y-direction 
%                    is imposed to a single horizontal layer of nodes
%                    if == 0, then no deformation is imposed on the nodes
% flag_nofail_zone : if == 1, then a no-fail zone is introduced by defining
%                    a mask (Boolean) array with a value of 1 for cells
%                    for which bonds connected to them are prevented to break
%                    (cells in the no-fail zone) and 0 for the other cells.
%                    Recall each cell has a node associated to it which shares 
%                    the same number. The cells in the no-fail zone are those 
%                    with associated nodes from the single horizontal layer of 
%                    nodes being deformed (if flag_deformation == 1) having
%                    positive x-coordinates in the underformed configuration
%                    if == 0, a no-fail zone is not introduced          

% Discussion:
% ----------
% We create a uniform grid within the domain (-1,1)x(-1,1) and plot
% the grid in the current configuration (which can be either undeformed or 
% deformed) with all the intact bonds, based on the FA algorithm.
%
% We define a set of candidate nodes on which to impose displacements: the 
% first horizontal layer of nodes above the x-axis. These nodes are highlighted
% (in red) in the plot. 
%
% For the case where we impose a deformation, we break all bonds that exceed 
% a critical stretch, unless they are connected to a point in a no-fail zone
% (if a no-fail zone is introduced).
%
% Setting the critical stretch: 
% 
% The condition for bond breaking is (see BondBreaking function): 
% 
%   rk2  > ((so+1)*Rk)^2. 
% 
% Define the critical current bond length as 1.25*h with h a uniform grid 
% spacing, and assume all bonds have reference length h (i.e., m-ratio = 1). 
% Then, the corresponding critical stretch is obtained by setting rk2 =
% (1.25*h)^2 and Rk = h, and replacing ">" by "=" in the bond-breaking condition: 
%   
%   (1.25*h)^2 = ((so+1)*h)^2 => 1.25*h = (so+1)*h => so = 1.25 - 1 = 0.25.

function PlotBondBreaking(flag_deformation,flag_nofail_zone)

    % Current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                        Material constants 
    % --------------------------------------------------------------------

    % Horizon
    del = 0.2;

    % Critical stretch
    so = 0.25;

    % Influence function order indicator
    omega = 0;

    % Algorithm for computation of neighbor areas
    AlgName = 'FA';

    % --------------------------------------------------------------------
    %                     Domain and discretization
    % --------------------------------------------------------------------

    % Domain boundaries
    Xo = -1;   Xn = 1;
    Yo = -1;   Yn = 1;

    % Grid perturbation coefficient
    PG = 0;

    % Uniform grid spacing
    h = del;
    dx = h;
    dy = dx;

    % Number of nodes in the x-direction
    Nx = round((Xn - Xo)/dx);
    % Number of nodes in the y-direction
    Ny = round((Yn - Yo)/dy);

    % --------------------------------------------------------------------
    %                          Generate grid
    % --------------------------------------------------------------------

    [xx,yy,~,~,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);

    % --------------------------------------------------------------------
    %                        Create neighbor list
    % --------------------------------------------------------------------

    % Flag for Rectangular Domain Uniform Grid (RDUG)
    flag_RDUG = 1;

    [u_NA,~,~,r_hat_NA,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);

    % --------------------------------------------------------------------
    %                        Impose deformation 
    % --------------------------------------------------------------------

    % Initialize displacements
    v = 0*xx;
    w = 0*yy;

    % Define set of nodes to impose displacements: 1st horizontal layer
    % of nodes above x-axis
    I_disp = find(yy > 0 & yy < h);

    if flag_deformation == 1
        % Impose deformation
        w(I_disp) = (0.5*h)*(1-abs(xx(I_disp)));
    elseif flag_deformation == 0

    else
        % Change directory
        cd(Directory)

        error('flag_deformation should be 0 or 1.')
    end

    % --------------------------------------------------------------------
    %               Define no-fail mask (Boolean) array:
    % 1st horizontal layer of nodes with positive x-coordinates above x-axis 
    % --------------------------------------------------------------------

    if flag_nofail_zone == 1
        % No-fail function
        nofailfunc = @(x,y) ( y > 0 & y < h & x > 0 ); 
        % No-fail mask (Boolean) array
        mask_nofail = nofailfunc(xx,yy);
    elseif flag_nofail_zone == 0
        % No-fail mask (Boolean) array
        mask_nofail = 0*xx;
    else
        % Change directory
        cd(Directory)

        error('flag_nofail_zone should be 0 or 1.')
    end

    % --------------------------------------------------------------------
    %                           Break bonds
    % --------------------------------------------------------------------

    [u_NA] = BondBreaking(xx,yy,v,w,so,u_NA,r_hat_NA,x_hat_NA,y_hat_NA,mask_nofail);

    % --------------------------------------------------------------------
    %                            Plot grid 
    %                     (current configuration)
    % --------------------------------------------------------------------

    % Change directory
    cd(Directory)

    % Plot grid
    PlotGrid(xx+v,yy+w,M,xx1,yy1,'b','w')
    box('on')

    % --------------------------------------------------------------------
    %             Highlight nodes with imposed displacements
    % --------------------------------------------------------------------

    hold all

    scatter(xx(I_disp)+v(I_disp),yy(I_disp)+w(I_disp),'r','filled')

    % --------------------------------------------------------------------
    %                     Plot all intact bonds 
    % --------------------------------------------------------------------

    % Find maximum number of neighbors any node can have
    zmax = length(u_NA(1,:));

    % Number of nodes
    Nnodes = length(xx);

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
                line([xi + v(ui) xk + v(uk)],[yi + w(ui) yk + w(uk)],'Color','c','LineStyle','-')

            end 
        end
    end

    box('on')

end
