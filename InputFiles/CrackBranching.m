
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
% Input deck for a crack branching problem on a pre-notched soda–lime glass 
% thin plate 
% ========================================================================

% Reference: 
% ---------
% F. Bobaru and G. Zhang, Why do cracks branch? A peridynamic investigation of
% dynamic brittle fracture, International Journal of Fracture 196 (2015): 59-98.

% ------------------------------------------------------------------------
%                  Domain geometry and discretization
% ------------------------------------------------------------------------

% Domain boundaries
Xo = -0.05; % [m] : Left  boundary of the domain 
Xn =  0.05; % [m] : Right boundary of the domain
Yo = -0.02; % [m] : Lower boundary of the domain
Yn =  0.02; % [m] : Upper boundary of the domain

% Number of nodes
Nx = 300; % Number of nodes in the x-direction 
Ny = 120; % Number of nodes in the y-direction

% Grid perturbation coefficient
PG = 0;

% ------------------------------------------------------------------------
%                         Time discretization
% ------------------------------------------------------------------------

% Initial time
Ti = 0;

% Final time
Tf = 4.3e-5; % [s] :    43 microsec

% Time step   
dt = 6.7E-8; % [s] : 67E-3 microsec

% Time-integration scheme
TimeScheme = 'VVerlet'; % Velocity Verlet 

% ------------------------------------------------------------------------
%                              PD model 
% ------------------------------------------------------------------------

% Constitutive model
model = 'GPMB'; % Generalized Prototype Microelastic Brittle (GPMB)

% Plane elasticity model
PlanarModel = 'PlaneStress';

% Horizon
del = 0.001; % [m]

% Influence function order indicator
omega = 0; % Constant influence function

% Flag for bond breaking
flag_BB = 1;

% ------------------------------------------------------------------------
%                       Classical material properties
% ------------------------------------------------------------------------

% Mass density 
rho = 2440; % [kg/m^3]

% Young's modulus   
E = 72e+9;  % [Pa] : 72 GPa 

% Fracture energy 
Go = 3.8;   % [J/m^2]

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

% Compute dy
dy = (Yn - Yo)/Ny;

% Traction amplitude
sigma = 2E6; % [Pa] : 2 MPa

bvfunc = @(x,y,t) (0.*x + 0.*y)*t;                         % x-component of body force density
bwfunc = @(x,y,t) ( abs(y) > Yn - dy ).*sigma.*sign(y)/dy; % y-component of body force density

% ------------------
% Initial conditions
% ------------------ 

% Initial displacement functions

vofunc = @(x,y) (0.*x + 0.*y); % x-component of initial displacement
wofunc = @(x,y) (0.*x + 0.*y); % y-component of initial displacement

% Initial velocity functions

Vvofunc = @(x,y) (0.*x + 0.*y); % x-component of initial velocity
Vwofunc = @(x,y) (0.*x + 0.*y); % y-component of initial velocity

% -------------------
% No-fail zone
% -------------------

% No-fail function
nofailfunc = @(x,y) ( abs(y) > Yn - del ); 

% ------------------
% Prenotch
% ------------------

%                       Xc1   Yc1  Xc2  Yc2
PreNotchCoordinates = [-0.05  0.0  0.0  0.0];

% ------------------------------------------------------------------------
%                           Postprocessing
% ------------------------------------------------------------------------

% Flag for plotting during time integration
flag_DynamicPlotting = 1;

% Frequency of plotting during time integration
DynamicPlotFrequency = 40; % Plot every 40 time steps (beginning from the first one)

% Frequency of time-integration step display
TimeStepDisplayFrequency = 20;

% Flag for plotting at final time
flag_FinalPlots = 1;

% Plot settings
%                     Field Name      Field variable   Colorbar title  Point size  Colormap limits  Colormap    Axes limits  Configuration
PlotSettings = {'StrainEnergyDensity' , 'log10(W)'  , '$\log_{10}(W)$'  ,    8    ,     [0 3.5]   ,   'jet'  , [Xo Xn Yo Yn] , 'Reference';
                      'Damage'        ,   'phi'     ,   '$\varphi$'     ,    8    ,     [0 0.4]   , 'parula' , [Xo Xn Yo Yn] , 'Reference'};

% Flag to visualize pre-notch as damaged
flag_DamagedPrenotches = 1;

% Flag to create video(s): works only if flag_DynamicPlotting = 1
flag_video = 1;

% Video frame rate
video_frate = 20;

% ------------------------------------------------------------------------
