
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
% The function NeighborAreaBondLengthCoord computes the area associated 
% with a neighboring quadrature point, the corresponding bond length, and 
% the coordinates of the quadrature point
% ========================================================================

% Input
% -----
% xi      : x coordinate of node i (source node)
% xk      : x coordinate of node k (neighbor node)
% yi      : y coordinate of node i (source node)
% yk      : y coordinate of node k (neighbor node)
% dx      : mesh spacing in the x-direction
% dy      : mesh spacing in the y-direction
% xx1k    : x coordinates of the vertices of cell k
% yy1k    : y coordinates of the vertices of cell k
% del     : horizon
% VVk     : area in reference configuration of cell k
% AlgName : name of algorithm (in string format) for computation of neighbor areas
%           'FA'    : FA algorithm
%           'PA-AC' : PA-AC algorithm
%           'IPA-AC': IPA-AC algorithm

% Output
% ------
% Vk      : area associated with neighbor cell k relative to node i 
% rk_hat  : bond length in reference configuration between node i and
%           the quadrature point associated with the cell k relative to node i
% xk_hat  : x-coordinate of the quadrature point associated with the cell k relative to node i
% yk_hat  : y-coordinate of the quadrature point associated with the cell k relative to node i

% Discussion
% ----------
% Vk is the quadrature weight: 
% - if AlgName = 'FA', Vk is simply the full area (FA) of the cell k; 
% - if AlgName = 'PA-AC' or 'IPA-AC', Vk is the area of the overlapping 
%   region between the neighborhood of node i and the cell k (assumed to be
%   a square cell), which can be a partial area (PA).
% 
% rk_hat is the effective bond length: 
% - if AlgName = 'FA' or 'PA-AC', rk_hat is the "standard" bond length, i.e., 
%   the distance between node i and node k;
% - if AlgName = 'IPA-AC', rk_hat is the distance between node i and the centroid 
%   of the overlapping region between the neighborhood of node i and the cell k
%   (assumed to be a square cell).
%
% (xk_hat,yk_hat) are the coordinates of the quadrature point:
% - if AlgName = 'FA' or 'PA-AC', the quadrature point is the node k; 
% - if AlgName = 'IPA-AC', the quadrature point is the centroid of the overlapping 
%   region between the neighborhood of node i and the cell k (assumed to be a square cell).
%
% A partial area (PA) is computed exactly following an analytical
% calculation (AC) in the reference below (see Algorithm PA-AC in page 192). 
% An improved (I) algorithm that can compute both a PA and its associated 
% centroid is also presented in the reference below (see Algorithm IPA-AC in page 198). 

% Reference: 
% ---------
% P. Seleson, Improved one-point quadrature algorithms for two-dimensional
% peridynamic models based on analytical calculations, Computer Methods in 
% Applied Mechanics and Engineering 282 (2014): 184–217.

function [Vk,rk_hat,xk_hat,yk_hat] = NeighborAreaBondLengthCoord(xi,xk,yi,yk,dx,dy,xx1k,yy1k,del,VVk,AlgName)

    % --------------------------------------------------------------------
    %                    Check validity of inputs
    % --------------------------------------------------------------------

    if dx <= 0  
        error('dx should be positive.')
    end

    if dy <= 0  
        error('dy should be positive.')
    end

    if length(xx1k) ~= length(yy1k) 
        error('xx1k and yy1k have different lengths.')
    end

    if del <= 0  
        error('del should be positive.')
    end

    % --------------------------------------------------------------------

    % Tolerance
    tol = 1E-15;

    % Compute distance squared between node i and node k
    r2 = ((xk-xi)^2+(yk-yi)^2);

    % Check that node k is not the same as node i
    if r2 < tol
        error('Node k should be different from node i.')
    end
    
    if strcmp(AlgName,'FA')

        % ----------------------------------------------------------------
        %                        FA algorithm
        % ----------------------------------------------------------------

        if r2 < del^2 + tol
            % Neighbor area: full area
            Vk = VVk;
            % Bond length
            rk_hat  = sqrt(r2);
        else
            % Neighbor area
            Vk = 0;
            % Bond length
            rk_hat  = 0;
        end

        % Assign the quadrature point coordinates to be those of node k
        xk_hat = xk;
        yk_hat = yk;

    elseif strcmp(AlgName,'PA-AC') || strcmp(AlgName,'IPA-AC')

        % ----------------------------------------------------------------
        %              PA-AC algorithm or IPA-AC algorithm 
        % ----------------------------------------------------------------

        % Check if grid is uniform
        if abs(dx-dy) > tol 
            error('PA-AC and IPA-AC algorithms require a uniform grid.')
        end

        % Define uniform grid spacing
        h = dx;

        % Assign cell k vertices
        xvertex = xx1k;
        yvertex = yy1k;

        % ------------------------------------
        % Check if cell k is a square of side  
        % h around node k 
        % ------------------------------------

        % Check if cell k has four vertices
        if length(xvertex) ~= 4 
            error('PA-AC and IPA-AC algorithms require a square cell k.')
        end

        % Check if the four cell k vertices form a square cell of side h around node k
        minxk = xk - 0.5*h;
        maxxk = xk + 0.5*h;
        minyk = yk - 0.5*h;
        maxyk = yk + 0.5*h;

        if abs(xvertex(1) - minxk) > tol || abs(xvertex(2) - maxxk) > tol || abs(xvertex(3) - maxxk) > tol || abs(xvertex(4) - minxk) > tol
            error('Cell k is not a square cell of side dx around node k.')
        end

        if abs(yvertex(1) - maxyk) > tol || abs(yvertex(2) - maxyk) > tol || abs(yvertex(3) - minyk) > tol || abs(yvertex(4) - minyk) > tol
            error('Cell k is not a square cell of side dx around node k.')
        end

        % ------------------------------------
        
        % Initialize flags for cell mapping to first quadrant
        flag_lr     = 0; % Left-right mirroring
        flag_bt     = 0; % Bottom-top mirroring
        flag_rotate = 0; % Rotation

        % ------------------------------------
        % Map all points to first quadrant
        % ------------------------------------
        
        % Mirror cell k left to right
        if xk < xi - tol
            xk = xi + (xi- xk);
            xvertex(1) = xk - 0.5*h; 
            xvertex(2) = xvertex(1) + h;
            xvertex(3) = xvertex(2);
            xvertex(4) = xvertex(1);
            flag_lr = 1;
        end

        % Mirror cell k bottom to top
        if yk < yi - tol
            yk = yi + (yi - yk);
            yvertex(1) = yk + 0.5*h; 
            yvertex(2) = yvertex(1);
            yvertex(3) = yvertex(1) - h;
            yvertex(4) = yvertex(3);
            flag_bt = 1;
        end

        % Rotate cell k from +y-axis to +x-axis (clockwise rotation)
        if abs(xk - xi) < tol 
            deltaY = yk - yi;
            yk = yi;
            xk = xi + deltaY;

            xvertex(1) = xk - 0.5*h; 
            xvertex(2) = xvertex(1) + h;
            xvertex(3) = xvertex(2);
            xvertex(4) = xvertex(1);
            yvertex(1) = yk + 0.5*h; 
            yvertex(2) = yvertex(1);
            yvertex(3) = yvertex(1) - h;
            yvertex(4) = yvertex(3); 
            flag_rotate = 1;
        end
        
        % Initialize the quadrature point coordinates to be those of 
        % (possibly transformed) node k
        xk_hat = xk;
        yk_hat = yk;

        % ------------------------------------
        % Count number of corners of cell k 
        % inside neighborhood of node i
        % ------------------------------------

        counter = 0;
        for i = 1:4
            if (xvertex(i) - xi)^2 + (yvertex(i) - yi)^2 < del^2 + tol
                counter = counter + 1;
            end
        end

        % ------------------------------------
        % Compute neighbor area and associated
        % centroid (for IPA-AC Algorithm)
        % ------------------------------------ 

        if counter == 4

            % -----------------------------------------------
            %                    Case I
            % -----------------------------------------------
            % Neighbor area: full area
            Vk = VVk;
            % Bond length
            rk_hat = sqrt(r2);
            % -----------------------------------------------

        elseif counter == 3

            % -----------------------------------------------
            %                    Case II
            % ----------------------------------------------- 
            
            H1 = yk + 0.5*h - yi;
            L2 = xk + 0.5*h - xi;

            H2 = sqrt(del^2 - L2^2);
            L1 = sqrt(del^2 - H1^2);

            d = sqrt( (H1 - H2)^2 + (L2 - L1)^2 );
            l = sqrt( del^2 - (0.5*d)^2 );
            gamma = asin(0.5*d / del);

            % Neighbor area: partial area
            A1 = h*(h - (H1 - H2));
            A2 = (h - (L2 - L1))*(H1 - H2);
            A3 = 0.5*(L2 - L1)*(H1 - H2);
            A4 = gamma*(del^2) - 0.5*d*l;

            Vk = A1 + A2 + A3 + A4;

            % Bond length
            if strcmp(AlgName,'PA-AC')
                rk_hat = sqrt(r2);
            else
                % Centroid
                x1 = xi + L2 - 0.5*h;
                y1 = yi + H2 - 0.5*(h - (H1 - H2));
                x2 = xi + L1 - 0.5*(h - (L2 - L1));
                y2 = yi + H1 - 0.5*(H1 - H2);
                x3 = xi + L2 - (2/3)*(L2 - L1);
                y3 = yi + H1 - (2/3)*(H1 - H2);

                theta = atan(H2/L2);
                l_hat = (4*del*(sin(gamma)^3)) / (3*(2*gamma - sin(2*gamma)));
                x4 = xi + l_hat*cos(theta + gamma);
                y4 = yi + l_hat*sin(theta + gamma);

                xk_hat = (A1*x1 + A2*x2 + A3*x3 + A4*x4) / Vk;
                yk_hat = (A1*y1 + A2*y2 + A3*y3 + A4*y4) / Vk;

                rk_hat = sqrt((xk_hat-xi)^2+(yk_hat-yi)^2);
            end
            % -----------------------------------------------

        elseif counter == 2

            if abs(yk - yi) < tol

                if xk + 0.5*h > xi + del - tol

                    % -----------------------------------------------
                    %                  Case III(a2)
                    % -----------------------------------------------  

                    l = sqrt(del^2 - (0.5*h)^2);
                    gamma = asin(0.5*h/del);

                    % Neighbor area: partial area
                    A1 = (l - (xk - 0.5*h - xi))*h;
                    A2 = gamma*del^2 - 0.5*h*l;

                    Vk = A1 + A2;

                    % Bond length
                    if strcmp(AlgName,'PA-AC')    
                        rk_hat = sqrt(r2);
                    else
                        % Centroid
                        x1 = xi + l - 0.5*(l - (xk - 0.5*h - xi));
                        l_hat  = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );
                        x2 = xi + l_hat;

                        xk_hat = (A1*x1 + A2*x2)/Vk;
                        yk_hat = yi;

                        rk_hat = xk_hat-xi;
                    end
                    % -----------------------------------------------

                else

                    % -----------------------------------------------
                    %                  Case III(b)
                    % -----------------------------------------------  

                    l = sqrt(del^2 - (0.5*h)^2);
                    L = xk + 0.5*h - xi;

                    gamma = acos(l/del);

                    if abs(del-L) > tol
                        beta  = acos(L/del);
                        d = 2*sqrt(del^2 - L^2);
                    else
                        beta  = acos(1);
                        d = 0;
                    end

                    % Neighbor area: partial area
                    A1 = (h - (L-l))*h;
                    A2 = gamma*del^2 - 0.5*h*l;
                    A3 = beta*del^2  - 0.5*d*L;

                    Vk = A1 + A2 - A3;

                    % Bond length
                    if strcmp(AlgName,'PA-AC')    
                        rk_hat = sqrt(r2);
                    else
                        % Centroid
                        x1 = xi + l - 0.5*(h - (L-l));
                        l_hat  = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );
                        x2 = xi + l_hat;
                        l_hat_prime  = (4*del*(sin(beta))^3) / ( 3*(2*beta - sin(2*beta)) );
                        x3 = xi + l_hat_prime;

                        xk_hat = (A1*x1 + A2*x2 - A3*x3)/Vk;
                        yk_hat = yi;

                        rk_hat = xk_hat-xi;
                    end
                    % -----------------------------------------------

                end

            else

                rbr2 = (xvertex(3) - xi)^2 + (yvertex(3) - yi)^2;

                if rbr2 > del^2 

                    % -----------------------------------------------
                    %                  Case III(a1)
                    % -----------------------------------------------  

                    H1 = yk + 0.5*h - yi;
                    H2 = yk - 0.5*h - yi;

                    L1 = sqrt(del^2 - H1^2);
                    L2 = sqrt(del^2 - H2^2);

                    d = sqrt( (L2-L1)^2 + h^2 );
                    l = sqrt(del^2 - (0.5*d)^2);
                    gamma = asin(0.5*d/del);

                    % Neighbor area: partial area
                    A1 = (L1 - (xk - 0.5*h - xi))*h;
                    A2 = 0.5*(L2-L1)*h;
                    A3 = gamma*del^2 - 0.5*d*l;

                    Vk = A1 + A2 + A3;

                    % Bond length
                    if strcmp(AlgName,'PA-AC')    
                        rk_hat = sqrt(r2);
                    else
                        % Centroid
                        x1 = xi + L1 - 0.5*(L1 - (xk - 0.5*h - xi));
                        y1 = yi + H1 - 0.5*h;
                        x2 = xi + L2 - 2*(L2-L1)/3;
                        y2 = yi + H1 - 2*h/3;
                        theta = atan(H2/L2);
                        l_hat = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );
                        x3 = xi + l_hat*cos(theta+gamma);
                        y3 = yi + l_hat*sin(theta+gamma);

                        xk_hat = (A1*x1 + A2*x2 + A3*x3)/Vk;
                        yk_hat = (A1*y1 + A2*y2 + A3*y3)/Vk;

                        rk_hat = sqrt((xk_hat-xi)^2+(yk_hat-yi)^2);
                    end
                    % -----------------------------------------------

                else

                    % -----------------------------------------------
                    %                  Case III(c)
                    % -----------------------------------------------  

                    L1 = xk - 0.5*h - xi;
                    L2 = xk + 0.5*h - xi;

                    H1 = sqrt(del^2 - L1^2);
                    H2 = sqrt(del^2 - L2^2);

                    d = sqrt(h^2 + (H1-H2)^2);
                    l = sqrt(del^2 - (0.5*d)^2);
                    gamma = asin(0.5*d/del);

                    % Neighbor area: partial area
                    A1 = (H2 - (yk-0.5*h-yi))*h;
                    A2 = 0.5*h*(H1-H2);
                    A3 = gamma*del^2 - 0.5*d*l;

                    Vk = A1 + A2 + A3;

                    % Bond length
                    if strcmp(AlgName,'PA-AC')    
                        rk_hat = sqrt(r2);
                    else
                        % Centroid
                        x1 = xi + L2 - 0.5*h;
                        y1 = yi + H2 - 0.5*(H2 - (yk-0.5*h-yi));
                        x2 = xi + L1 + h/3;
                        y2 = yi + H2 + (H1-H2)/3;
                        theta = atan(H2/L2);
                        l_hat  = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );
                        x3 = xi + l_hat*cos(theta+gamma);
                        y3 = yi + l_hat*sin(theta+gamma);

                        xk_hat = (A1*x1 + A2*x2 + A3*x3)/Vk;
                        yk_hat = (A1*y1 + A2*y2 + A3*y3)/Vk;

                        rk_hat = sqrt((xk_hat-xi)^2+(yk_hat-yi)^2);
                    end
                    % -----------------------------------------------

                end

            end

        elseif counter == 1

            % -----------------------------------------------
            %                   Case IV
            % ----------------------------------------------- 

            L1 = xk - 0.5*h - xi;
            H2 = yk - 0.5*h - yi;

            H1 = sqrt(del^2 - L1^2);
            L2 = sqrt(del^2 - H2^2);

            d = sqrt( (L2-L1)^2 + (H1-H2)^2);
            l = sqrt(del^2 - (0.5*d)^2);
            gamma = asin(0.5*d/del);

            % Neighbor area: partial area
            A1 = 0.5*(L2-L1)*(H1-H2);
            A2 = gamma*del^2 - 0.5*d*l;

            Vk = A1 + A2;

            % Bond length
            if strcmp(AlgName,'PA-AC')    
                rk_hat = sqrt(r2);
            else
                % Centroid
                x1 = xi + L1 + (L2-L1)/3;
                y1 = yi + H2 + (H1-H2)/3;
                theta = atan(H2/L2);
                l_hat  = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );
                x2 = xi + l_hat*cos(theta+gamma);
                y2 = yi + l_hat*sin(theta+gamma);

                xk_hat = (A1*x1 + A2*x2)/Vk;
                yk_hat = (A1*y1 + A2*y2)/Vk;

                rk_hat = sqrt((xk_hat-xi)^2+(yk_hat-yi)^2);
            end
            % -----------------------------------------------

        else

            if abs(yk - yi) < tol && xk - 0.5*h < xi + del

                % -----------------------------------------------
                %                     Case V
                % -----------------------------------------------   

                l = xk - 0.5*h - xi;

                if abs(del-l) > tol
                    d = 2*sqrt(del^2 - l^2);
                    gamma = acos(l/del);
                else
                    d = 0;
                    gamma = acos(1);
                end

                % Neighbor area: partial area
                Vk = gamma*del^2 - 0.5*d*l;

                % Bond length
                if strcmp(AlgName,'PA-AC')    
                    rk_hat = sqrt(r2);
                else
                    % Centroid
                    if abs(del-l) > tol
                        l_hat  = (4*del*(sin(gamma))^3) / ( 3*(2*gamma - sin(2*gamma)) );

                        xk_hat = xi + l_hat;
                        yk_hat = yi;

                        rk_hat = xk_hat-xi;
                    else
                        xk_hat = xi + delta;
                        yk_hat = yi;

                        rk_hat = xk_hat-xi;
                    end
                end
                % -----------------------------------------------

            else
                Vk = 0;
                rk_hat  = 0;
            end

        end
        
        % Move centroid back to original quarter
        if( flag_rotate == 1)
            deltaX = xk_hat - xi;
            xk_hat = xi;
            yk_hat = yi + deltaX;
        end
        
        if( flag_bt == 1)
            yk_hat = 2*yi - yk_hat;
        end
        
        if( flag_lr == 1)
            xk_hat = 2*xi - xk_hat;
        end

    else

        error('Invalid AlgName.')

    end

end
