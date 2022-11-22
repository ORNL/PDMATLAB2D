
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
% The function GridGenerator generates the computational grid for simulation
% ========================================================================

% Input
% -----
% Xo  : left boundary of the domain
% Xn  : right boundary of the domain
% Yo  : lower boundary of the domain
% Yn  : upper boundary of the domain
% Nx  : number of nodes in the x-direction 
% Ny  : number of nodes in the y-direction
% PG  : grid perturbation coefficient - PG in [0,1) - 0 is for no perturbation
%       (perturbation of up to 0.25*dx in the x-direction 
%        and up to 0.25*dy in the y-direction)

% Output
% ------
% xx  : x coordinates of all nodes in the grid (1D array of length Nx x Ny)
% yy  : y coordinates of all nodes in the grid (1D array of length Nx x Ny)
% x   : array of coordinates in the x-direction (1D array of length Nx)
% y   : array of coordinates in the y-direction (1D array of length Ny)
% dx  : mesh spacing in the x-direction
% dy  : mesh spacing in the y-direction
% VV  : area of each cell corresponding to each node (1D array of length Nx x Ny)
% xx1 : x coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% yy1 : y coordinates of all cell vertices (1D array of length (Nx+1) x (Ny+1))
% M   : array that maps each node to its cell vertices 
%       (2D array of size (Nx x Ny)-by-4)

function [xx,yy,x,y,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG)

    % --------------------------------------------------------------------
    %                    Check validity of inputs
    % --------------------------------------------------------------------

    if Xn <= Xo  
        error('Xn should be larger than Xo.')
    end

    if Yn <= Yo
        error('Yn should be larger than Yo.')
    end

    if floor(Nx) ~= ceil(Nx) || Nx <= 0
        error('Nx is not a natural number.')
    end

    if floor(Ny) ~= ceil(Ny) || Ny <= 0
        error('Ny is not a natural number.')
    end

    if PG < 0 || PG >= 1
        error('Invalid PG.')
    end
    
    % --------------------------------------------------------------------
    %                  Create x,y and x1,y1 vectors
    % --------------------------------------------------------------------
    %
    % Note: We count nodes in the following order: left to right, top to bottom. 
    %       Consequently, we flip the y vector.

    % Compute dx and dy
    dx = (Xn - Xo)/Nx;
    dy = (Yn - Yo)/Ny;

    % x and y vectors for inner domain (cell center nodes)
    x = Xo+dx/2:dx:Xn-dx/2;
    y = Yo+dy/2:dy:Yn-dy/2;
    y_flipped = fliplr(y);

    % x and y vectors for full domain (cell vertices)
    x1 = Xo:dx:Xn;
    y1 = Yo:dy:Yn;
    y1_flipped = fliplr(y1);
    
    Nx1 = length(x1);
    Ny1 = length(y1);
    
    % --------------------------------------------------------------------
    %          Create connectivity matrix M and all-nodes 
    %               and all-vertices x and y vectors
    % --------------------------------------------------------------------
    % Construct M matrix such that each row represents a node's 4 vertices
    % (i.e., the vertices of the cell represented by that node).
    % The 4 entries in each row correspond to the top-left, top-right
    % bottom-right, and bottom-left vertices, respectively:
    %
    %  col1       col2
    %    ----------
    %   |          |
    %   |          |
    %   |          |
    %   |          |
    %    ----------
    %  col4       col3
    %
    %   The grid nodes are counted as follows:
    % 
    %    .1      .2     .3   ...  .Nx
    %   
    %    .(Nx+1) .(Nx+2)    ...   .(2*Nx)
    %    .                          .
    %    .                          .
    %    .                          .
    %    .(Ny-1)*Nx+1     ...     .Nx*Ny
      
    % Initialize meshgrid corresponding to nodes and vertices

    % Nodes:
    xx = zeros(Nx*Ny,1);
    yy = zeros(Nx*Ny,1);

    % Vertices:
    xx1 = zeros(Nx1*Ny1,1);
    yy1 = zeros(Nx1*Ny1,1);

    % Initialize nodal areas
    VV = zeros(Nx*Ny,1);
    
    % Initialize columns of M matrix that maps each node to its 4 vertices
    col1 = zeros(Nx*Ny,1);
    col2 = zeros(Nx*Ny,1);
    col3 = zeros(Nx*Ny,1);
    col4 = zeros(Nx*Ny,1);
    
    % Loop over cells
    count = 0;
    for ii = 1:Ny   % Loop top to bottom (nodes)
        for jj = 1:Nx  % Loop left to right (nodes)
            count = count + 1;
            col1(count) = (ii-1)*Nx1 + jj;
            col2(count) = (ii-1)*Nx1 + jj + 1;
            col4(count) = ii*Nx1 + jj;
            col3(count) = ii*Nx1 + jj + 1;
        end
    end
    
    M = [col1,col2,col3,col4];
    [s1,~] = size(M);
    
    % Loop through vertices and fill vertex vectors xx1 and yy1
    for i = 1:Ny1    % Loop top to bottom (vertices)
        for j = 1:Nx1   % Loop left to right (vertices)         
          	% Convert (i,j) into vertex counter number
         	u1 = (i-1)*Nx1 + j;
            xx1(u1) = x1(j);
            yy1(u1) = y1_flipped(i);
        end
    end
    
    % --------------------------------------------------------------------
    %         Perturb vertices in the interior of the domain
    %            (leave vertices on boundary unperturbed)
    % --------------------------------------------------------------------
    %
    % Note: We perturb nodes by at most 0.25 times the grid spacing in each
    %       dimension to guarantee the cell remains convex
    
    if PG > 0
        for i = 2:Ny1-1    % Loop top to bottom (inner vertices)
            for j = 2:Nx1-1   % Loop left to right (inner vertices)             
                % Convert (i,j) into vertex counter number
                u1 = (i-1)*Nx1 + j;
                xx1(u1) = xx1(u1) + PG*0.25*dx*(-1+2*rand(1));
                yy1(u1) = yy1(u1) + PG*0.25*dy*(-1+2*rand(1));
            end
        end
    end
    
    % --------------------------------------------------------------------
    %            Compute areas and centroids of cells
    % --------------------------------------------------------------------

    if PG > 0
        % Loop through rows of M and calculate the area of each perturbed cell
        % as well as its centroid. Use these values to fill VV, xx, and yy.
        for u = 1:s1 % Loop over rows of M (same as number of cells)
            % Decompose each quadrilateral cell into two triangles:
            % Triangle 1: A-B1-C (col1-col2-col3)
            % Triangle 2: A-B2-C (col1-col4-col3)
            cVec = M(u,:);
            A   = cVec(1);
            B1  = cVec(2);
            C   = cVec(3);
            B2  = cVec(4);
            
            % Compute areas of cells
            AT1 = abs((xx1(A)*(yy1(B1)-yy1(C))+xx1(B1)*(yy1(C)-yy1(A))+xx1(C)*(yy1(A)-yy1(B1)))/2);
            AT2 = abs((xx1(A)*(yy1(B2)-yy1(C))+xx1(B2)*(yy1(C)-yy1(A))+xx1(C)*(yy1(A)-yy1(B2)))/2);
            VV(u) = AT1 + AT2;
            
            % Compute centroids of cells
            Ox1 = (xx1(A)+xx1(B1)+xx1(C)) / 3;
            Ox2 = (xx1(A)+xx1(B2)+xx1(C)) / 3;
            Oy1 = (yy1(A)+yy1(B1)+yy1(C)) / 3;
            Oy2 = (yy1(A)+yy1(B2)+yy1(C)) / 3;
            xx(u) = (Ox1*AT1 + Ox2*AT2) / VV(u);
            yy(u) = (Oy1*AT1 + Oy2*AT2) / VV(u);
        end       
    else
        % Compute areas and centroids for regular grid
        
        % Compute areas of cells
        VV(:) = dx*dy;
        
        % Compute centroids of cells
        for i = 1:Ny   % Loop top to bottom (nodes)
            for j = 1:Nx   % Loop left to right (nodes)
                % Convert (i,j) into node number
                u = (i-1)*Nx + j;
                xx(u) = x(j);
                yy(u) = y_flipped(i);
            end
        end  
    end
    
end
