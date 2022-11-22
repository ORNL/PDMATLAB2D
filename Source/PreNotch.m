
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
% The function PreNotch creates a pre-notch described by an arbitrary line 
% segment 
% ========================================================================

% Input
% -----
% xx        : x coordinates of all nodes in the grid 
% yy        : y coordinates of all nodes in the grid 
% u_NA      : array of neighbor numbers for all nodes
% (Xc1,Yc1) : coordinates of one endpoint of the pre-notch ("crack")
% (Xc2,Yc2) : coordinates of the other endpoint of the pre-notch ("crack")

% Output
% ------
% u_NA      : updated array of neighbor numbers for all nodes for the creation 
%             of the pre-notch

% Discussion:
% ----------
% A pre-notch is created by defining a line segment representing it, and then
% breaking all bonds that either cross the line segment or overlap it.
%
% The function employs the function SegmentsIntersection to checks whether 
% two segments intersect

function [u_NA] = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2)
    
    % Tolerance
    tol = 1E-15;
   
    % --------------------------------------------------------------------
    %                    Check validity of inputs
    % --------------------------------------------------------------------

    % Check if xx and yy arrays have the same length
    if length(xx) ~= length(yy)
        error('Arrays xx and yy have different lengths.')
    end

    % Number of nodes
    Nnodes = length(xx);

    % Size of u_NA
    [s1,~] = size(u_NA);

    % Check consistency of u_NA array dimensions with number of nodes
    if s1 ~= Nnodes 
        error('u_NA is inconsistent with the number of nodes.')
    end

    % Check if pre-notch has finite length
    if abs(Xc1 - Xc2) < tol && abs(Yc1-Yc2) < tol
        error('Pre-notch has zero length.')
    end

    % --------------------------------------------------------------------
    %              Break bonds that intersect pre-notch
    % --------------------------------------------------------------------

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

            % Get neighbor cell number uk
            uk = u_NA(ui,z);
    
            % Only consider neighbor cells
            if uk > 0

                % Get x-coordinate of node uk
                xk = xx(uk);
    
                % Get y-coordinate of node uk
                yk = yy(uk);
    
                % Check if bond ui-uk intersects pre-notch
                flag_intersection = SegmentsIntersection(xi,yi,xk,yk,Xc1,Yc1,Xc2,Yc2);
    
                if flag_intersection == 1

                    % Remove cell uk from neighbor list of node ui
                    u_NA(ui,z) = 0;

                    % Read neighbor list of node uk
                    uVeck = u_NA(uk,:);

                    % Find where cell ui is in the neighbor list of node uk
                    mask = (uVeck==ui);
                    
                    % Remove cell ui from neighbor list of node uk
                    u_NA(uk,mask) = 0;

                end
            end
        end
    end
end
