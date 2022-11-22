
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
% The function PlotSegmentsIntersection plots two random line segments,  
% showing one of the two segments in different colors depending on whether 
% an intersection occurs
% ========================================================================

% Input
% -----
% ncase : case number  
%         1: Collinear line segments
%         2: Parallel but not collinear line segments
%         3: Non-parallel line segments 
%         4: Non-parallel lines segments with 2nd line segment intersecting 
%            1st line segment endpoint

function PlotSegmentsIntersection(ncase)

    % Identify current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                         Define domain
    % --------------------------------------------------------------------

    % Domain limits in the x-direction
    Xo = -0.5;
    Xn =  0.5;

    % Domain limits in the y-direction
    Yo = -0.5;
    Yn =  0.5;
    
    % --------------------------------------------------------------------
    %                      1st random line segment
    % --------------------------------------------------------------------
    
    % Coordinates of the one end of the 1st line segment 
    Xp1A = Xo + rand(1,1)*(Xn-Xo); 
    Yp1A = Yo + rand(1,1)*(Yn-Yo); 

    % Coordinates of the other end of the 1st line segment  
    Xp2A = Xo + rand(1,1)*(Xn-Xo); 
    Yp2A = Yo + rand(1,1)*(Yn-Yo); 

    % Define line containing 1st segment: (vector) p + t * (vector) r 
    p = [Xp1A Yp1A];
    r = [Xp2A-Xp1A,Yp2A-Yp1A];

    % --------------------------------------------------------------------
    %                      2nd random line segment
    % --------------------------------------------------------------------

    switch ncase

        case 1

            % ------------------------------------------------------------
            %             Case 1: Collinear line segments
            % ------------------------------------------------------------

            % Define endpoints of 2nd line segment: q1 & q2
            u  = -0.5 + rand(1,2);
            q1 = p  + u(1)*r;
            q2 = q1 + u(2)*r;

            % Coordinates of the one end of the 2nd line segment 
            Xp1B = q1(1);
            Yp1B = q1(2);

            % Coordinates of the other end of the 2nd line segment  
            Xp2B = q2(1);
            Yp2B = q2(2);

        case 2

            % ------------------------------------------------------------
            %       Case 2: Parallel but not collinear line segments
            % ------------------------------------------------------------

            % Define endpoints of 2nd line segment: q1 & q2
            u  = -0.5 + rand(1,2);
            q1 = p  + u(1)*r;
            q2 = q1 + u(2)*r;

            % Shift q1 and q2 by a perpendicular vector to r
            r_perp = 0.5*rand(1,1)*[-r(2),r(1)]./norm(r);
            q1 = q1 + r_perp;
            q2 = q2 + r_perp;

            % Coordinates of the one end of the 2nd line segment 
            Xp1B = q1(1);
            Yp1B = q1(2);

            % Coordinates of the other end of the 2nd line segment  
            Xp2B = q2(1);
            Yp2B = q2(2);

        case 3

            % ------------------------------------------------------------
            %             Case 3: Non-parallel line segments      
            % ------------------------------------------------------------
            
            % Randomly select 2nd line segment endpoints

            % Coordinates of the one end of the 2nd line segment
            Xp1B = -1 + 2*rand(1,1);
            Yp1B = -1 + 2*rand(1,1);

            % Coordinates of the other end of the 2nd line segment
            Xp2B = -1 + 2*rand(1,1);
            Yp2B = -1 + 2*rand(1,1);

        case 4

            % -------------------------------------------------------------
            %            Case 4: Non-parallel lines segments 
            % with 2nd line segment intersecting 1st line segment endpoint
            % -------------------------------------------------------------

            % Define unit vector in random orientation
            n = rand(1,2);
            n = n/norm(n);

            % Define endpoints of 2nd line segment: q1 & q2
            u = -0.5 + rand(1,2);
            q1 = p - abs(u(1))*n;
            q2 = p + abs(u(2))*n;

            Xp1B = q1(1);
            Yp1B = q1(2);
            Xp2B = q2(1);
            Yp2B = q2(2);

        otherwise

            % Change directory
            cd(Directory)

            error('Invalid case number.');

    end

    % --------------------------------------------------------------------
    %                 Check line segments intersection
    % --------------------------------------------------------------------

    flag_intersection = SegmentsIntersection(Xp1A,Yp1A,Xp2A,Yp2A,Xp1B,Yp1B,Xp2B,Yp2B);
  
    % --------------------------------------------------------------------
    %                       Plot line segments
    % --------------------------------------------------------------------

    % Note: 
    % - 1st line segment is plotted in black
    % - 2nd line segment is plotted in red if it intersects the 1st one
    %   and in blue otherwise

    % Create figure
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    
    % Plot 1st line segment
    plot([Xp1A Xp2A],[Yp1A Yp2A],'o-k','LineWidth',2)
     
     % Plot 2nd line segment
    if flag_intersection == 1
        plot([Xp1B Xp2B],[Yp1B Yp2B],'s-r','LineWidth',2);
    else
        plot([Xp1B Xp2B],[Yp1B Yp2B],'s-b','LineWidth',2);
    end

    % Define figure axes to guarantee 2nd line segment is within the domain
    axis([-1 1 -1 1])
   
    % Set axes font size and latex font style
    set(axes1,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex')

    % Axes labels
    xlabel('$x$','Interpreter','latex','FontSize',30)
    ylabel('$y$','Interpreter','latex','FontSize',30)

    % Box
    box('on') 

    % Change directory
    cd(Directory)

end
