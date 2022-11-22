
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
% The function InfluenceFunction computes influence functions
% ========================================================================

% Input
% -----
% omega : influence function order indicator
%         0  : constant = 1
%         0.5: piecewise constant
%         1  : piecewise linear
%         3  : piecewise cubic
%         5  : piecewise quintic
%         7  : piecewise septic
% r     : distance value(s)
% del   : horizon

% Output
% ------
% IF    : influence function value(s) for given distance value(s)

% Note
% ----
% The piecewise influence functions are such that they have a polynomial
% functional form for r <= del and are 0 otherwise.

function [IF] = InfluenceFunction(omega,r,del)

    % --------------------------------------------------------------------
    %                    Check validity of inputs
    % --------------------------------------------------------------------

    if ~all(r >= 0)
        error('r should contain nonnegative values only.');
    end

    if del <= 0  
        error('del should be positive.')
    end

    % --------------------------------------------------------------------

    % Tolerance
    tol = 1E-15;

    % Evaluate influence function (IF)
    if omega == 0
        % Constant IF
        IF = ones(size(r));
    else
        % Initialize IF
        IF = zeros(size(r));
     
        % Identify vector indexes for which r <= delta
        mask = (r < del + tol);
        
        switch omega
            case 0.5
                % Piecewise constant IF
                IF(mask) = 1;
            case 1
                % Piecewise linear IF
                IF(mask) = 1 - (r(mask)./del);
            case 3
                % Piecewise cubic IF
                IF(mask) = 1 - 3*(r(mask)./del).^2 + 2*(r(mask)./del).^3;
            case 5
                % Piecewise quintic IF
                IF(mask) = 1 - 10*(r(mask)./del).^3 + 15*(r(mask)./del).^4 - 6*(r(mask)./del).^5;
            case 7
                % Piecewise septic IF
                IF(mask) = 1 - 35*(r(mask)./del).^4 + 84*(r(mask)./del).^5 - 70*(r(mask)./del).^6 + 20*(r(mask)./del).^7;
            otherwise
                error('Invalid omega.');
        end
    end

end
