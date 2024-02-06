
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
% The function ForceEnergyDensity computes the internal force density and 
% the macroelastic energy density for all nodes
% ========================================================================

% Input
% -----
% xx         : x coordinates of all nodes in the grid 
% yy         : y coordinates of all nodes in the grid 
% v          : displacement of each node in the x-direction
% w          : displacement of each node in the y-direction
% c          : micromodulus constant
% u_NA       : array of neighbor numbers for all nodes
% IF_NA      : array of influence function values of neighbor bonds for all nodes 
% V_NA       : array of neighbor areas for all nodes 
% r_hat_NA   : array of reference lengths of neighbor bonds for all nodes 
% x_hat_NA   : array of x-coordinates of quadrature points for all nodes
% y_hat_NA   : array of y-coordinates of quadrature points for all nodes 
% model      : constitutive model 'GPMB'
%              'GPMB' : generalized prototype microelastic brittle
% flag_RDUG  : if == 1, then the grid is uniform over a rectangular domain
%              (dx = dy); the flag name RDUG stands for "Rectangular Domain
%              Uniform Grid" 
%              if == 0, then the grid is a general grid

% Output
% ------
% Fv         : x-component of internal force density for all nodes
% Fw         : y-component of internal force density for all nodes
% W          : macroelastic energy density for all nodes

% Discussion:
% ----------
% The GPMB model is a generalization of the PMB model by incorporating an 
% influence function. 
% 
% The PMB model was presented in:
%
% S.A. Silling and E. Askari, A meshfree method based on the peridynamic 
% model of solid mechanics, Computers and Structures 83 (2005): 1526–1535.
% 
% The GPMB model was presented in:
%
% P. Seleson and M. L. Parks, On the role of the influence function in the 
% peridynamic theory, International Journal for Multiscale Computational 
% Engineering 9(6) (2011): 689–706. 
%
% The input variable "model" currently only takes 'GPMB' as valid input. 
% Additional models can be added using the overall if statement by adding
% an elseif statement to specify a condition for each desired model. As an example, 
% if the 'X' model is needed, one could extend the if statement as follows:
%
% if strcmp(model,'GPMB')
%
% elseif strcmp(model,'X')
%
% else
%
%     error('Invalid model.')
%
% end

function [Fv,Fw,W] = ForceEnergyDensity(xx,yy,v,w,c,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,model,flag_RDUG)

    % --------------------------------------------------------------------
    %                           GPMB model
    % --------------------------------------------------------------------

    if strcmp(model,'GPMB')
    
        % Number of nodes
        Nnodes = length(xx);
    
        % Initialize internal force density components and macroelastic energy density arrays
        Fv = zeros(Nnodes,1); % Array of x-components of internal force density for all nodes
        Fw = zeros(Nnodes,1); % Array of y-components of internal force density for all nodes
        W  = zeros(Nnodes,1); % Array of macroelastic energy density for all nodes
   
        % Find maximum number of neighbors any node can have
        zmax = length(u_NA(1,:));
    
        if flag_RDUG == 1

            % --------------------------------------------------------
            %              Computation for uniform grid
            %                 over rectangular domain
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells with a larger number
                    % than the source node to compute the bond pairwise force
                    % and pairwise potential once per bond
                    if uk > ui

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = c*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % Update internal force function components of node uk
                        Fv(uk) = Fv(uk) - fvk*Vk;
                        Fw(uk) = Fw(uk) - fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*c*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                        % Update macroelastic energy density of node uk
                        W(uk) = W(uk) + 0.5*wk*Vk;

                    end
                end
            end

        elseif flag_RDUG == 0    

            % --------------------------------------------------------
            %              Computation for general grids
            % --------------------------------------------------------

            % Loop through all nodes ui
            for ui = 1:Nnodes

                % Get x-coordinate of node ui
                xi = xx(ui);
                % Get y-coordinate of node ui
                yi = yy(ui);

                % Get x-component of displacement of node ui
                vi = v(ui);
                % Get y-component of displacement of node ui
                wi = w(ui);

                % Loop over neighbors of node ui
                for z = 1:zmax

                    % Get neighbor cell number
                    uk = u_NA(ui,z);

                    % Only consider neighbor cells
                    if uk > 0

                        % Get x-coordinate of quadrature point
                        xk_hat = x_hat_NA(ui,z);
                        % Get y-coordinate of quadrature point
                        yk_hat = y_hat_NA(ui,z);

                        % Get x-component of displacement of node uk
                        vk = v(uk);
                        % Get y-component of displacement of node uk
                        wk = w(uk);

                        % Read per-bond quantities for bond ui-uk
                        IFk = IF_NA(ui,z);         % Influence function
                        Vk  = V_NA(ui,z);          % Neighbor area
                        Rk  = r_hat_NA(ui,z);      % Reference bond length

                        % Compute deformed bond components and length
                        rxk = xk_hat-xi + vk-vi;   % x-component of current relative position
                        ryk = yk_hat-yi + wk-wi;   % y-component of current relative position
                        rk  = sqrt(rxk^2 + ryk^2); % current bond length

                        % Compute bond stretch
                        sk = (rk - Rk)/Rk;

                        % ------------------------------------------------
                        %        Compute internal force density
                        % ------------------------------------------------

                        % Compute bond pairwise force magnitude
                        fk = c*IFk*sk;

                        % Compute bond pairwise force components
                        fvk = fk*rxk/rk;
                        fwk = fk*ryk/rk;

                        % Update internal force function components of node ui
                        Fv(ui) = Fv(ui) + fvk*Vk;
                        Fw(ui) = Fw(ui) + fwk*Vk;

                        % ------------------------------------------------
                        %      Compute macroelastic energy density
                        % ------------------------------------------------

                        % Compute bond pairwise potential
                        wk = 0.5*c*IFk*(sk^2)*Rk;

                        % Update macroelastic energy density of node ui
                        W(ui) = W(ui) + 0.5*wk*Vk;

                    end
                end
            end

        else

            error('flag_RDUG should be 0 or 1.')

        end

    else

        error('Invalid model.')

    end

end
