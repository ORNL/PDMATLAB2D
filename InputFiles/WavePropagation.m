
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
% Input deck for a wave propagation problem with an initial displacement 
% based on a radial Gaussian distribution (the parameters are provided in
% consistent units)
% ========================================================================

% ------------------------------------------------------------------------
%                  Domain geometry and discretization
% ------------------------------------------------------------------------

% Domain boundaries
Xo = 0; % Left  boundary of the domain
Xn = 5; % Right boundary of the domain
Yo = 0; % Lower boundary of the domain
Yn = 5; % Upper boundary of the domain

% Number of nodes
Nx = 200; % Number of nodes in the x-direction 
Ny = 200; % Number of nodes in the y-direction

% Grid perturbation coefficient
PG = 0;

% ------------------------------------------------------------------------
%                         Time discretization
% ------------------------------------------------------------------------

% Initial time
Ti = 0;

% Final time
Tf = 1.8;

% Time step
dt = 0.01;

% Time-integration scheme
TimeScheme = 'VVerlet'; % Velocity Verlet

% ------------------------------------------------------------------------
%                             PD model 
% ------------------------------------------------------------------------

% Constitutive model
model = 'GPMB'; % Generalized Prototype Microelastic Brittle (GPMB) model

% Plane elasticity model
PlanarModel = 'PlaneStress';

% Horizon
del = 0.1; 

% Influence function order indicator
omega = 0; % Constant influence function

% Flag for bond breaking
flag_BB = 0;

% ------------------------------------------------------------------------
%                      Classical material properties
% ------------------------------------------------------------------------

% Mass density
rho = 1;

% Young's modulus
E = 1;

% Fracture energy
Go = 1;

% ------------------------------------------------------------------------
%                         Meshfree discretization 
% ------------------------------------------------------------------------

% Algorithm for computation of neighbor areas
AlgName = 'FA'; % FA algorithm

% ------------------------------------------------------------------------
%                           Problem settings
% ------------------------------------------------------------------------

% ------------------
% Body force density
% ------------------

bvfunc = @(x,y,t) (0.*x + 0.*y)*t; % x-component of body force density
bwfunc = @(x,y,t) (0.*x + 0.*y)*t; % y-component of body force density

% ------------------
% Initial conditions
% ------------------ 

% Initial displacement functions

% Tolerance
tol = 1E-15;

% Parameters
xm    = 2.5;     % x-coordinate of pulse center 
ym    = 2.5;     % y-coordinate of pulse center 
A     = 0.025;   % amplitude of radial Gaussian distribution
sigma = 1/30;    % standard deviation of radial Gaussian distribution
mu    = 6*sigma; % radial distance from pulse center of radial Gaussian distribution mean 

% Functions
r      = @(x,y) sqrt((x-xm).^2+(y-ym).^2);                    % distance from pulse center
uo     = @(x,y) ( r(x,y) > mu-6*sigma & r(x,y) < mu+6*sigma ).* A .* exp( (-(r(x,y) - mu).^2) / (2*sigma)^2 ); % magnitude of initial displacement
vofunc = @(x,y) uo(x,y).*(x-xm)./(r(x,y) + tol);              % x-component of initial displacement
wofunc = @(x,y) uo(x,y).*(y-ym)./(r(x,y) + tol);              % y-component of initial displacement

% Initial velocity functions

Vvofunc = @(x,y) 0.*x + 0.*y; % x-component of initial velocity
Vwofunc = @(x,y) 0.*x + 0.*y; % y-component of initial velocity

% ------------------------------------------------------------------------
%                           Postprocessing
% ------------------------------------------------------------------------

% Flag for plotting during time integration
flag_DynamicPlotting = 1;

% Frequency of plotting during time integration
DynamicPlotFrequency = 10; % Plot every 10 time steps (beginning from the first one)

% Frequency of time-integration step display
TimeStepDisplayFrequency = 10;

% Flag for plotting at final time
flag_FinalPlots = 1;

% Plot settings
%                     Field Name             Field variable         Colorbar title   Point size  Colormap limits   Colormap    Axes limits    Configuration
PlotSettings = {'DisplacementMagnitude' , 'sqrt(v.^2 + w.^2)/A'  , '$\|{\bf u}\|/A$' ,    8    ,     [0 1.0]   ,   'parula'  , [Xo Xn Yo Yn] , 'Reference'};

% Flag to create video(s): works only if flag_DynamicPlotting = 1
flag_video = 1;

% Video frame rate
video_frate = 20;

% ------------------------------------------------------------------------