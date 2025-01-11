
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
% Main script for running a PDMATLAB2D simulation
% ========================================================================

% Check if simulation output directory exists and create it otherwise
if ~exist(['../Outputs/' InputDeck], 'dir')
    mkdir(['../Outputs/' InputDeck])
end

% ------------------------------------------------------------------------
%                       Create video file(s)
% ------------------------------------------------------------------------

if flag_DynamicPlotting == 1

    % Check if video flag is defined
    if exist('flag_video','var')

        if flag_video == 1
            % Check if video output directory exists and create it otherwise
            if ~exist(['../Videos/' InputDeck], 'dir')
                mkdir(['../Videos/' InputDeck])
            end

            % Find number of plots: a video is generated for each plot
            [s1,~] = size(PlotSettings);

            % Initialize VideoWriter array
            vidfile = VideoWriter.empty(s1, 0);

            % Loop over plots
            for nplot = 1:s1
                % Read plot field name
                field_name = PlotSettings{nplot,1};

                % Video file name
                video_filename = ['../Videos/' InputDeck '/' field_name '.mp4'];

                % Open video file
                vidfile(nplot) = VideoWriter(video_filename,'MPEG-4');
                vidfile(nplot).FrameRate = video_frate;
                open(vidfile(nplot));
            end
        end

    else

        % Define null video flag
        flag_video = 0;

    end

end

% ------------------------------------------------------------------------
%                    Generate grid and neighbor list
% ------------------------------------------------------------------------

% Check if GridFile variable to load grid is defined
if exist('GridFile','var')

    % Check if grid file exists
    if exist(GridFile,'file') == 2

        % ---------------------------------
        %    Load grid and neighbor list
        % ---------------------------------
        tic
        load(GridFile);
        fprintf('Load grid and neighbor list ..................... = %f (sec) \n',toc)

        % Flag for Rectangular Domain Uniform Grid (RDUG)
        flag_RDUG = 0;

        % ---------------------------------
        %        Check grid inputs
        % ---------------------------------
        % Check if variables exist
        gridvars = ["xx","yy","u_NA","IF_NA","V_NA","r_hat_NA","x_hat_NA","y_hat_NA"];

        for n = 1:length(gridvars)           
            if exist(gridvars(n),'var') == 0
                cd ..
                error('Loaded grid variable %s does not exist.', gridvars(n))
            end
        end

        % Check consistency of xx and yy array dimensions
        if isequal(size(xx),size(yy))

        else
            cd ..
            error('xx and yy arrays have inconsistent dimensions.')
        end

        % Check consistency of u_NA, IF_NA, V_NA, r_hat_NA, x_hat_NA, and y_hat_NA array dimensions
        if isequal(size(u_NA),size(IF_NA),size(V_NA),size(r_hat_NA),size(x_hat_NA),size(y_hat_NA)) 

        else  
            cd ..
            error('u_NA, IF_NA, V_NA, r_hat_NA, x_hat_NA, and y_hat_NA arrays do not have all consistent dimensions.')        
        end

    else

        cd ..
        error('Grid file does not exist.')    
    
    end

else
    
    % ---------------------------------
    %        Generate grid
    % ---------------------------------

    tic
    [xx,yy,x,y,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);
    fprintf('Generate grid ................................... = %f (sec) \n',toc)

    % Tolerance
    tol = 1E-15;

    % Flag for Rectangular Domain Uniform Grid (RDUG)
    if abs(dx - dy) < tol 
        flag_RDUG = 1;
    else
        flag_RDUG = 0;
    end

    % ---------------------------------
    %     Generate neighbor list
    % ---------------------------------

    tic
    [u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);
    fprintf('Generate neighbor list .......................... = %f (sec)\n',toc)

end

% ------------------------------------------------------------------------
%                        Create no-fail mask
% ------------------------------------------------------------------------

if exist('nofailfunc','var')
    mask_nofail = nofailfunc(xx,yy);
else
    mask_nofail = 0*xx;
end

% ------------------------------------------------------------------------
%                       Create pre-notch(es)
% ------------------------------------------------------------------------

% Check if array of pre-notches coordinates is defined
if exist('PreNotchCoordinates','var')

    tic

    % Find number of pre-notches
    [s1,~] = size(PreNotchCoordinates);

    % Loop over pre-notches
    for n = 1:s1
        % Coordinates of one endpoint of the pre-notch
        Xc1 = PreNotchCoordinates(n,1);
        Yc1 = PreNotchCoordinates(n,2);

        % Coordinates of the other endpoint of the pre-notch
        Xc2 = PreNotchCoordinates(n,3);
        Yc2 = PreNotchCoordinates(n,4);

        % Create pre-notch
        [u_NA] = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);
    end

    fprintf('Create pre-notch(es) ............................ = %f (sec)\n',toc)

end

% ------------------------------------------------------------------------
%                   Compute peridynamic constants
% ------------------------------------------------------------------------

tic

[c,so] = PDBondConstants(omega,del,E,Go,model,PlanarModel);

fprintf('Compute PD constants ............................ = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                      Impose initial conditions
% ------------------------------------------------------------------------

tic

% Compute initial displacement for all nodes
v = vofunc(xx,yy);   % x-component of initial displacement
w = wofunc(xx,yy);   % y-component of initial displacement

% Compute initial velocity for all nodes
Vv = Vvofunc(xx,yy); % x-component of initial velocity
Vw = Vwofunc(xx,yy); % y-component of initial velocity

fprintf('Impose initial conditions ....................... = %f (sec)\n',toc)

% ------------------------------------------------------------------------
% Compute initial internal force density and macroelastic energy density
% ------------------------------------------------------------------------

tic

[Fv,Fw,W] = ForceEnergyDensity(xx,yy,v,w,c,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,model,flag_RDUG);

fprintf('Compute initial internal force & energy densities = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                Compute initial body force density
% ------------------------------------------------------------------------

tic

bv = bvfunc(xx,yy,0); % x-component of initial body force density
bw = bwfunc(xx,yy,0); % y-component of initial body force density

fprintf('Compute initial body force density .............. = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                 Compute denominator of damage ratio
% ------------------------------------------------------------------------

% Check if bond breaking is enabled
if flag_BB == 0

elseif flag_BB == 1
    
    if exist('flag_DamagedPrenotches','var')

        % Compute damage ratio denominator before application of pre-notches
        if flag_DamagedPrenotches == 1

            phiD = sum(V_NA,2);
        
        % Compute damage ratio denominator after application of pre-notches
        elseif flag_DamagedPrenotches == 0
        
            phiD = sum((u_NA>0).*V_NA,2);
        
        else
        
            error('flag_DamagedPrenotches should be 0 or 1.')
        
        end

    else

        % Compute default damage ratio denominator (after application of pre-notches)
        phiD = sum((u_NA>0).*V_NA,2);
    
    end

else

    error('flag_BB should be 0 or 1.')

end

% ------------------------------------------------------------------------
%                         Time integration loop
% ------------------------------------------------------------------------

% Create time vector
tVec = Ti:dt:Tf;

% Number of time steps
Nt = length(tVec);

% Initial time
t = Ti;

tic

% Loop over time steps
for n = 1:Nt-1

    [v,w,Vv,Vw,Fv,Fw,bv,bw,W,u_NA] = TimeIntegrator(TimeScheme,xx,yy,v,w,Vv,Vw,Fv,Fw,bv,bw,t,bvfunc,bwfunc,dt,u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA,rho,c,model,flag_RDUG,so,mask_nofail,flag_BB);

    % Time-integration step display
    if mod(n,TimeStepDisplayFrequency) == 0
        fprintf('Time integration (n = %4g / %g ) = %f (sec)\n',n,Nt-1,toc)
    end

    % Time after performing time-integration step
    t = tVec(n+1);

    % ---------------------------------
    %         Plot field(s)
    % ---------------------------------

    if flag_DynamicPlotting == 0

    elseif flag_DynamicPlotting == 1

        if n == 1 || mod(n-1,DynamicPlotFrequency) == 0

            % Find number of plots
            [s1,~] = size(PlotSettings);

            % Loop over plots
            for nplot = 1:s1

                % Figure number
                nfig = nplot;                       
                
                % Read plot field name
                cnodes_name = PlotSettings{nplot,1};    

                % Check if compute damage for plotting
                if flag_BB == 0

                elseif flag_BB == 1
                    if strcmp(cnodes_name,'Damage')
                        % Compute damage
                        phi = 1 - sum((u_NA>0).*V_NA,2)./phiD;
                    end
                else
                    error('flag_BB should be 0 or 1.')
                end

                % Read plot field variable
                cnodes = eval(PlotSettings{nplot,2}); 

                % Read plot settings
                ctitle        = PlotSettings{nplot,3};       % Colorbar title
                psize         = PlotSettings{nplot,4};       % Point size
                climits       = PlotSettings{nplot,5};       % Colormap limits
                cmap          = PlotSettings{nplot,6};       % Colormap
                box           = PlotSettings{nplot,7};       % Axes limits
                configuration = PlotSettings{nplot,8};       % Configuration: 'Reference' or 'Current'

                % Plot field
                if strcmp(configuration,'Reference')
                    PlotField(nfig,xx,yy,cnodes,ctitle,psize,climits,cmap,box)
                elseif strcmp(configuration,'Current')
                    PlotField(nfig,xx+v,yy+w,cnodes,ctitle,psize,climits,cmap,box)
                else
                    error('Invalid configuration.')
                end

                % Video update
                if flag_video == 1
                    % Access figure window
                    f = figure(nfig);

                    % Format time for title
                    time = sprintf('%.2e', t);
                    time = regexprep(time, '(e[\+\-])0(\d)', '$1$2');

                    % Set figure title
                    title(['Time = ' time '\,s'], 'Interpreter', 'latex')

                    % Capture frame and write it to video file
                    frame = getframe(gcf);
                    writeVideo(vidfile(nplot),frame);
                end

            end

        end

    else

        error('flag_DynamicPlotting should be 0 or 1.')

    end

end

% ------------------------------------------------------------------------
%                           Final outputs
% ------------------------------------------------------------------------

if flag_FinalPlots == 0

elseif flag_FinalPlots == 1
    
    % Find number of plots
    [s1,~] = size(PlotSettings);

    % Loop over plots
    for nplot = 1:s1

        % Figure number
        nfig = nplot;

        % Read plot field name
        cnodes_name = PlotSettings{nplot,1};

        % Check if compute damage for plotting
        if flag_BB == 0

        elseif flag_BB == 1
            if strcmp(cnodes_name,'Damage')
                % Compute damage
                phi = 1 - sum((u_NA>0).*V_NA,2)./phiD;
            end
        else
            error('flag_BB should be 0 or 1.')
        end

        % Read plot field variable
        cnodes = eval(PlotSettings{nplot,2});

        % Read plot settings
        ctitle        = PlotSettings{nplot,3};       % Colorbar title
        psize         = PlotSettings{nplot,4};       % Point size
        climits       = PlotSettings{nplot,5};       % Colormap limits
        cmap          = PlotSettings{nplot,6};       % Colormap
        box           = PlotSettings{nplot,7};       % Axes limits
        configuration = PlotSettings{nplot,8};       % Configuration: 'Reference' or 'Current'

        % Plot field
        if strcmp(configuration,'Reference')
            PlotField(nfig,xx,yy,cnodes,ctitle,psize,climits,cmap,box)
        elseif strcmp(configuration,'Current')
            PlotField(nfig,xx+v,yy+w,cnodes,ctitle,psize,climits,cmap,box)
        else
            error('Invalid configuration.')
        end

        % Save figure
        filename = ['../Outputs/' InputDeck '/' cnodes_name '.eps'];
        saveas(gcf,filename,'epsc');

        if flag_DynamicPlotting == 1

            % Video update
            if flag_video == 1
                % Access figure window
                f = figure(nfig);

                % Format time for title
                time = sprintf('%.2e', t);
                time = regexprep(time, '(e[\+\-])0(\d)', '$1$2');

                % Set figure title
                title(['Time = ' time '\,s'], 'Interpreter', 'latex')

                % Capture frame and write it to video file
                frame = getframe(gcf);
                writeVideo(vidfile(nplot),frame);

                % Close video file
                close(vidfile(nplot))
            end
            
        end

    end

else

    error('flag_FinalPlots should be 0 or 1.')

end
  