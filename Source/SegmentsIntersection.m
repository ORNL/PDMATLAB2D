
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
% The function SegmentsIntersection checks whether two line segments intersect
% ========================================================================

% Input
% -----
% (Xp1A,Yp1A) : coordinates of the  one  end of the 1st line segment (segment "A")
% (Xp2A,Yp2A) : coordinates of the other end of the 1st line segment (segment "A")
% (Xp1B,Yp1B) : coordinates of the  one  end of the 2nd line segment (segment "B")
% (Xp2B,Yp2B) : coordinates of the other end of the 2nd line segment (segment "B")

% Output
% ------
% flag_intersection: if == 1, then the line segments intersect
%                    if == 0, then the line segments do not intersect

% Discussion
% ----------
%
% The calculations in this function are based on an algorithm presented in:
% https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
% (accessed on November 6, 2022), which is a 2D adaptation of the reference below. 
% 
% The notation for the two lines is as follows:
%
% The 1st line (containing segment A) is defined by: (vector) p + t * (vector) r 
% The 2nd line (containing segment B) is defined by: (vector) q + u * (vector) s
%
% The function employs the function vcross2D (defined below) to compute 
% the "cross product" of two 2D vectors

% Reference: 
% ---------
% Goldman, R.: Intersection of two lines in three-space. In: Glassner, A.
% (ed.) Graphics Gems, p. 304. Academic Press, Inc., San Diego, CA (1990).
% Chap. V.3

function flag_intersection = SegmentsIntersection(Xp1A,Yp1A,Xp2A,Yp2A,Xp1B,Yp1B,Xp2B,Yp2B)

    % Tolerance
    tol = 1E-15;
    
    % Initialize flag for segments intersection
    flag_intersection = 0;

    % --------------------------------------------------------------------
    %                          Define lines
    % --------------------------------------------------------------------
    
    % 1st line:
    p = [Xp1A Yp1A];
    r = [Xp2A-Xp1A,Yp2A-Yp1A];
    
    % 2nd line:
    q = [Xp1B Yp1B];
    s = [Xp2B-Xp1B,Yp2B-Yp1B];
    
    % Check if r = [0 0]
    if norm(r) < tol
        error('Invalid 1st line segment.');
    end
    
    % Check if s = [0 0]
    if norm(s) < tol
        error('Invalid 2nd line segment.');
    end

    % --------------------------------------------------------------------
    %                 Check if line segments intersect
    % --------------------------------------------------------------------
    
    % "Cross product" r x s
    cross_r_s = vcross2D(r,s);
    
    % Check if lines are parallel
    if abs(cross_r_s) < tol
        
        % Check if lines are collinear
        if abs(vcross2D((q - p),r)) < tol
            
            % The lines are collinear: express endpoints of 1st line segment 
            %                          in terms of the 2nd line equation
            
            % Denominator
            dot_s_s = dot(s,s);
            
            % Parameter of endpoint p (of 1st line segment)  : find uo s.t. p = q + uo*s    
            uo = dot((p - q),s)/dot_s_s;
            
            % Parameter of endpoint p+r (of 1st line segment): find u1 s.t. p+r = q + u1*s
            %                           => r = q - p + u1*s = q - (q + uo*s) + u1*s = (u1-uo)*s 
            u1 = uo + dot(r,s)/dot_s_s;
 
            % Check if the line segments overlap
            if min(uo,u1) > 1 || max(uo,u1) < 0    
                
                % Line segments do not overlap
                
            else
  
                % Line segments overlap
                flag_intersection = 1;
                            
            end
            
        else
            
            % The lines are parallel but not collinear (line segments do not overlap)
            
        end
          
    else
        
        % The lines are not parallel: solve for t & u
        %
        % Note: there is one point s.t. (2nd line) q + u*s = p + t*r (1st line)
        % 
        %       Apply "cross product" with s:   (q + u*s) x s = (p + t*r) x s
        %                                     => q x s = p x s + t*(r x s)
        %                                     => (q - p) x s = t*(r x s)
        %                                     => t = ((q - p) x s) / (r x s)
        %
        %       Apply "cross product" with r:   (q + u*s) x r = (p + t*r) x r
        %                                     => q x r + u*(s x r) = p x r
        %                                     => u = (p - q) x r / (s x r) = (q - p) x r / (r x s)
        
        % Parameter t of 1st line
        t = vcross2D((q - p),s)/cross_r_s;

        % Parameter u of 2nd line
        u = vcross2D((q - p),r)/cross_r_s;
         
        % Check if segments intersect: intersection occurs if 0 <= u,t <= 1
        if -tol < u && u < 1 + tol &&  -tol < t && t < 1 + tol
            
            % Line segments intersect
            flag_intersection = 1;
            
        else
            
            % Line segments are not parallel but they do not intersect
            
        end
       
    end
    
end

% ========================================================================
% The function vcross2D computes the "cross product" of two 2D vectors
% ========================================================================

% Input
% -----
% v1     : first vector with coordinates  [v1(1) v1(2)]
% v2     : second vector with coordinates [v2(1) v2(2)]

% Output
% ------
% vcross : "cross product" v1 x v2

% Discussion
% ----------
% The 2D vectors v1 and v2 are cast as 3D vectors with a z-component equal 0. 
% Then, the "cross product" v1 x v2 is given by the z-component of the standard 
% 3D cross product.

function vcross = vcross2D(v1,v2)

    % "Cross product"
    vcross = v1(1)*v2(2) - v1(2)*v2(1);
    
end
