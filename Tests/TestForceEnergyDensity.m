
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
% The function TestForceEnergyDensity tests the function ForceEnergyDensity
% ========================================================================

% Discussion:
% ----------
% We impose a displacement given by either: 
%     (1) isotropic extension or 
%     (2) quadratic deformation
% for two plane elasticity models: 
%     (i)  plane strain
%     (ii) plane stress
% and compute the internal force density components and macroelastic energy 
% density for a point in the bulk of a system for different influence functions 
% and neighbor areas algorithms.
%
% The function outputs the values for each case together with associated 
% reference continuum values.

function TestForceEnergyDensity

    % Current directory
    Directory = pwd;

    % Change directory
    cd ../Source/

    % --------------------------------------------------------------------
    %                    Material constants and model 
    % --------------------------------------------------------------------

    % Horizon
    del = 0.2;

    % Young's modulus
    E = 1;

    % Fracture energy
    Go = 1;

    % Model
    model = 'GPMB';

    % --------------------------------------------------------------------
    %                     Domain and discretization
    % --------------------------------------------------------------------

    % Domain boundaries
    Xo = 0;   Xn = 1;
    Yo = 0;   Yn = 1;

    % Grid perturbation coefficient
    PG = 0;

    % m-ratio
    m = 6;

    % Grid spacing
    dx = del/m;
    dy = dx;

    % Number of nodes in the x-direction
    Nx = round((Xn - Xo)/dx);
    % Number of nodes in the y-direction
    Ny = round((Yn - Yo)/dy);

    % --------------------------------------------------------------------
    %                     Deformation coefficients
    % --------------------------------------------------------------------

    % Coefficient of isotropic extension 
    sbar = 0.1;

    % Coefficient of quadratic displacement
    V_11 = 0.01;

    % --------------------------------------------------------------------
    % Run over both deformation cases and for each use different options:
    % different plane elasticity models, influence functions, and neighbor 
    % areas algorithms
    % --------------------------------------------------------------------

    % Run over two deformation cases: isotropic extension (ncase = 1) and
    %                                 quadratic deformation (ncase = 2)
    for ncase = 1:2

        fprintf('\n')

        if ncase == 1
            fprintf('====================================================================================== \n')
            fprintf('                                 Isotropic extension \n')
            fprintf('====================================================================================== \n')
        else
            fprintf('======================================================================================  \n')
            fprintf('                                Quadratic deformation \n')
            fprintf('======================================================================================  \n')
        end

        % Run over plane elasticity models: plane strain and plane stress
        for PlanarModel = ["PlaneStrain" "PlaneStress"]

            if strcmp(PlanarModel,'PlaneStrain')
                fprintf('                ------------------------------------------------------- \n')
                fprintf('                                    Plane strain \n')
                fprintf('                ------------------------------------------------------- \n')
            else
                fprintf('                ------------------------------------------------------- \n')
                fprintf('                                   Plane stress \n')
                fprintf('                ------------------------------------------------------- \n')
            end
            
            % Print output labels
            fprintf('omega     AlgName            F_i                      F              W_i         W  \n')

            % Run over influence functions
            for omega = [0 0.5 1 3 5 7]

                % Run over neighbor areas algorithms
                for AlgName = ["FA" "PA-AC" "IPA-AC"]

                    % ----------------------------------------------------
                    %   Continuum expressions for internal force density
                    %          and macroelastic energy density
                    % ----------------------------------------------------

                    % Isotropic extension
                    if ncase == 1

                        % Internal force density
                        Fv_continuum = 0;
                        Fw_continuum = 0;

                        % Macroelastic energy density
                        if strcmp(PlanarModel,'PlaneStrain')
                            % Plane strain
                            W_continuum = @(x) (8/5)*E*sbar^2 + 0*x;
                        else
                            % Plane stress
                            W_continuum = @(x) (3/2)*E*sbar^2 + 0*x;
                        end

                    % Quadratic deformation
                    elseif ncase == 2

                        % Internal force density
                        if strcmp(PlanarModel,'PlaneStrain')
                            % Plane strain
                            Fv_continuum = (12/5)*V_11*E;
                            Fw_continuum = 0;
                        else
                            % Plane stress
                            Fv_continuum = (9/4)*V_11*E;
                            Fw_continuum = 0;
                        end

                        % Macroelastic energy density for each influence function (IF)
                        switch omega
                            case 0
                                % Constant IF
                                Wfactor = @(x) 24*x.^2 + 3*del^2;
                            case 0.5
                                % Piecewise constant IF
                                Wfactor = @(x) 24*x.^2 + 3*del^2;
                            case 1
                                % Piecewise linear IF
                                Wfactor = @(x) 24*x.^2 + 2*del^2;
                            case 3
                                % Piecewise cubic IF
                                Wfactor = @(x) 24*x.^2 + (45/28)*del^2;
                            case 5
                                % Piecewise quintic IF
                                Wfactor = @(x) 24*x.^2 + (7/5)*del^2;
                            case 7
                                % Piecewise septic IF
                                Wfactor = @(x) 24*x.^2 + (14/11)*del^2;
                        end

                        if strcmp(PlanarModel,'PlaneStrain')
                            % Plane strain
                            W_continuum = @(x) (E*V_11^2/10)*Wfactor(x);
                        else
                            % Plane stress
                            W_continuum = @(x) (3*E*V_11^2/32)*Wfactor(x);
                        end

                    else
                        error('Unknown deformation case.')
                    end

                    % ----------------------------------------------------
                    %                    Generate grid
                    % ----------------------------------------------------

                    [xx,yy,~,~,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);

                    % ----------------------------------------------------
                    %                 Create neighbor list
                    % ----------------------------------------------------

                    % Flag for Rectangular Domain Uniform Grid (RDUG)
                    flag_RDUG = 1;

                    [u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);

                    % ----------------------------------------------------
                    %                  Impose deformation
                    % ----------------------------------------------------

                    % Isotropic extension
                    if ncase == 1
                        
                        v = sbar*xx;
                        w = sbar*yy;

                    % Quadratic deformation
                    else
                        
                        v = V_11*xx.^2;
                        w = 0*yy;

                    end

                    % ----------------------------------------------------
                    %               Read micromodulus constant
                    % ----------------------------------------------------

                    [c,~] = PDBondConstants(omega,del,E,Go,model,PlanarModel);

                    % ----------------------------------------------------
                    %            Compute internal force density 
                    %            and macroelastic energy density
                    % ----------------------------------------------------

                    % Number of nodes
                    Nnodes = length(xx);

                    % Compute internal force density components and macroelastic
                    % energy density for all nodes
                    [Fv,Fw,W] = ForceEnergyDensity(xx,yy,v,w,c,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,model,flag_RDUG);

                    % ----------------------------------------------------
                    %         Print values for a node in the bulk
                    % ----------------------------------------------------

                    % Run over all nodes
                    for ui = 1:Nnodes

                        % Check if x-coordinate of node ui is farther than delta from domain boundary
                        if xx(ui) > Xo + del && xx(ui) < Xn - del

                            % Check if y-coordinate of node ui is farther than delta from domain boundary
                            if yy(ui) > Yo + del && yy(ui) < Yn - del

                                fprintf('  %3g      %6s   (%+5.2e,%+5.2e)   (%+5.2e,%+5.2e)  %5.2e   %5.2e \n',omega,AlgName,Fv(ui),Fw(ui),Fv_continuum,Fw_continuum,W(ui),W_continuum(xx(ui)))

                                break

                            end

                        end

                    end

                end

            end

        end

    end

    % Change directory
    cd(Directory)

end
