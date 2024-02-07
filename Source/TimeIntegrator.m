
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
% The function TimeIntegrator performs a time integration step
% ========================================================================

% Input
% -----
% TimeScheme  : time-integration scheme 'VVerlet'
%               'VVerlet': velocity Verlet
% xx          : x coordinates of all nodes in the grid 
% yy          : y coordinates of all nodes in the grid 
% v           : displacement of each node in the x-direction (at time t)
% w           : displacement of each node in the y-direction (at time t)
% Vv          : x-component of velocity for all nodes (at time t)
% Vw          : y-component of velocity for all nodes (at time t)
% Fv          : x-component of internal force density for all nodes (at time t)
% Fw          : y-component of internal force density for all nodes (at time t)
% bv          : x-component of body force density for all nodes (at time t)
% bw          : y-component of body force density for all nodes (at time t)
% t           : time
% bvfunc      : function for the x-component of body force density
% bwfunc      : function for the y-component of body force density
% dt          : time step
% u_NA        : array of neighbor numbers for all nodes (at time t)
% IF_NA       : array of influence function values of neighbor bonds for all nodes 
% V_NA        : array of neighbor areas for all nodes 
% r_hat_NA    : array of reference lengths of neighbor bonds for all nodes 
% x_hat_NA    : array of x-coordinates of quadrature points for all nodes
% y_hat_NA    : array of y-coordinates of quadrature points for all nodes 
% rho         : mass density
% c           : micromodulus constant
% model       : constitutive model 'GPMB'
%               'GPMB' : generalized prototype microelastic brittle
% flag_RDUG   : if == 1, then the grid is uniform over a rectangular domain
%               (dx = dy); the flag name RDUG stands for "Rectangular Domain
%               Uniform Grid" 
%               if == 0, then the grid is a general grid
% so          : critical stretch
% mask_nofail : Boolean array with a value of 1 for cells in a no-fail zone
%               and 0 for the other cells (recall each cell has a node 
%               associated to it which shares the same number)
% flag_BB     : if == 1, then the function allows bond breaking
%               if == 0, then the function does not allow bond breaking

% Output
% ------
% v           : displacement of each node in the x-direction (at time t+dt)
% w           : displacement of each node in the y-direction (at time t+dt)
% Vv          : x-component of velocity for all nodes (at time t+dt)
% Vw          : y-component of velocity for all nodes (at time t+dt)
% Fv          : x-component of internal force density for all nodes (at time t+dt)
% Fw          : y-component of internal force density for all nodes (at time t+dt)
% bv          : x-component of body force density for all nodes (at time t+dt)
% bw          : y-component of body force density for all nodes (at time t+dt)
% W           : macroelastic energy density for all nodes (at time t+dt)    
% u_NA        : updated array of neighbor numbers for all nodes after possible bond breaking (at time t+dt)

% Discussion:
% ----------
% The input variable "TimeScheme" currently only takes 'VVerlet' as valid input. 
% Additional time-integration schemes can be added using the overall if statement 
% by adding an elseif statement to specify a condition for each desired
% time-integration scheme. As an example, if the 'X' time-integration scheme 
% is needed, one could extend the if statement as follows:
%
% if strcmp(TimeScheme,'VVerlet')
%
% elseif strcmp(TimeScheme,'X')
%
% else
%
%     error('Invalid TimeScheme.')
%
% end

function [v,w,Vv,Vw,Fv,Fw,bv,bw,W,u_NA] = TimeIntegrator(TimeScheme,xx,yy,v,w,Vv,Vw,Fv,Fw,bv,bw,t,bvfunc,bwfunc,dt,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,rho,c,model,flag_RDUG,so,mask_nofail,flag_BB)

    % --------------------------------------------------------------------
    %              Velocity Verlet time-integration scheme
    % --------------------------------------------------------------------

    if strcmp(TimeScheme,'VVerlet')

        % The algorithm is ( r:position; v:velocity; a:acceleration ):
        % r(t+dt) = r(t) + dt*v(t) + (dt^2/2)*a(t)
        % v(t+dt) = v(t) + dt*( (a(t) + a(t+dt))/2 )
        %
        % This is equivalent to: 
        % Step 1: v(t + dt/2) = v(t) + (dt/2)*a(t)
        % Step 2: r(t + dt)   = r(t) + dt*v(t + dt/2)
        % Step 3: v(t + dt)   = v(t + dt/2) + (dt/2)*a(t+dt)
        %
        % Note: updating current positions is equivalent to updating 
        %       displacements, so the implementation below uses displacements

        % Step 1: Compute velocity at t+dt/2 for all nodes 
        V12v = Vv + ((0.5*dt.*(Fv + bv))./rho); % x-component of velocity
        V12w = Vw + ((0.5*dt.*(Fw + bw))./rho); % y-component of velocity
 
        % Step 2: Compute displacement at t+dt for all nodes  
        v = v + dt.*V12v; % x-component of displacement
        w = w + dt.*V12w; % y-component of displacement

        % Update neighbor list if bond breaking occurs
        if flag_BB == 0

        elseif flag_BB == 1
            [u_NA] = BondBreaking(xx,yy,v,w,so,u_NA,r_hat_NA,x_hat_NA,y_hat_NA,mask_nofail);
        else
            error('flag_BB should be 0 or 1.')
        end
        
        % Compute internal force density and macroelastic energy density at t+dt for all nodes
        [Fv,Fw,W] = ForceEnergyDensity(xx,yy,v,w,c,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,model,flag_RDUG);

        % Compute body force density at t+dt for all nodes 
        bv = bvfunc(xx,yy,t+dt); % x-component of body force density
        bw = bwfunc(xx,yy,t+dt); % y-component of body force density

        % Step 3: Compute velocity at t+dt for all nodes 
        Vv = V12v + ((0.5*dt.*(Fv + bv))./rho);
        Vw = V12w + ((0.5*dt.*(Fw + bw))./rho);

    else

        error('Invalid TimeScheme.')

    end
    
end
