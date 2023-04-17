
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
% The function PDBondConstants computes peridynamic constants based on
% material constants from the classical local theory
% ========================================================================

% Input
% -----
% omega       : influence function order indicator
% del         : horizon
% E           : Young's modulus
% Go          : fracture energy
% model       : constitutive model 'GPMB'
%               'GPMB' : generalized prototype microelastic brittle
% PlanarModel : plane elasticity model
%               'PlaneStrain' : plane strain assumption
%               'PlaneStress' : plane stress assumption

% Output
% ------
% c           : micromodulus constant
% so          : critical stretch

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

function [c,so] = PDBondConstants(omega,del,E,Go,model,PlanarModel)

    % --------------------------------------------------------------------
    %                    Check validity of inputs
    % --------------------------------------------------------------------

    if del <= 0  
        error('del should be positive.')
    end

    if E <= 0  
        error('E should be positive.')
    end

    if Go <= 0  
        error('Go should be positive.')
    end

    % --------------------------------------------------------------------
    %                           GPMB model
    % --------------------------------------------------------------------

    if strcmp(model,'GPMB')
    
        % ----------------------------------------------------------------
        %                         Plane strain
        % ----------------------------------------------------------------

        if strcmp(PlanarModel,'PlaneStrain')
 
              % Material constants for each influence function (IF)
              switch omega
                 case 0
                     % Constant IF
                     c  = (48*E)/(5*pi*del.^3);
                     so = sqrt((5*pi*Go)/(12*E*del));
                  case 0.5
                      % Piecewise constant IF
                      c  = (48*E)/(5*pi*del.^3);
                      so = sqrt((5*pi*Go)/(12*E*del));
                 case 1
                     % Piecewise linear IF
                     c  = (192*E)/(5*pi*del.^3);
                     so = sqrt((25*pi*Go)/(48*E*del));
                 case 3
                     % Piecewise cubic IF
                     c  = (48*E)/(pi*del.^3);
                     so = sqrt((7*pi*Go)/(12*E*del));
                 case 5
                     % Piecewise quintic IF
                     c  = (1344*E)/(25*pi*del.^3);
                     so = sqrt((5*pi*Go)/(8*E*del));
                 case 7
                     % Piecewise septic IF
                     c  = (288*E)/(5*pi*del.^3);
                     so = sqrt((55*pi*Go)/(84*E*del));
                 otherwise
                     error('Invalid omega.');
              end

        % ----------------------------------------------------------------
        %                         Plane stress
        % ----------------------------------------------------------------
    
        elseif strcmp(PlanarModel,'PlaneStress')
    
            % Material constants for each influence function (IF)
            switch omega
                case 0
                    % Constant IF
                    c  = (9*E)/(pi*del.^3);
                    so = sqrt((4*pi*Go)/(9*E*del));
                case 0.5
                    % Piecewise constant IF
                    c  = (9*E)/(pi*del.^3);
                    so = sqrt((4*pi*Go)/(9*E*del));
                case 1
                    % Piecewise linear IF
                    c  = (36*E)/(pi*del.^3);
                    so = sqrt((5*pi*Go)/(9*E*del));
                case 3
                    % Piecewise cubic IF
                    c  = (45*E)/(pi*del.^3);
                    so = sqrt((28*pi*Go)/(45*E*del));
                case 5
                    % Piecewise quintic IF
                    c  = (252*E)/(5*pi*del.^3);
                    so = sqrt((2*pi*Go)/(3*E*del));
                case 7
                    % Piecewise septic IF
                    c  = (54*E)/(pi*del.^3);
                    so = sqrt((44*pi*Go)/(63*E*del));
                otherwise
                    error('Invalid omega.');
            end
    
        else
    
            error('Invalid PlanarModel.')
    
        end
    
    else
    
        error('Invalid model.')
    
    end

end
