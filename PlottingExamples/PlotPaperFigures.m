
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
% The function PlotPaperFigures plots the figures in the paper describing
% the PDMATLAB2D code
% ========================================================================

% Input
% -----
% nfig: figure identifier in string format (e.g., '3(a)' for Figure 3(a))

function PlotPaperFigures(nfig)

    % Plotting functions directory
    Directory = pwd;

    % ====================================================================
    %                            Figure 3
    % ====================================================================

    if strcmp(nfig,'3(a)')
    
        % ---------------------------------------------------------------
        %                           Figure 3(a)
        % ---------------------------------------------------------------

        % Change directory
        cd ../Source/

        % Generate uniform grid
        [xx,yy,~,~,~,~,~,xx1,yy1,M] = GridGenerator(0,2,0,1,10,10,0);
    
        % Change directory
        cd(Directory)
    
        % Plot grid
        PlotGrid(xx,yy,M,xx1,yy1,'b','k')
    
    elseif strcmp(nfig,'3(b)')

        % ---------------------------------------------------------------
        %                           Figure 3(b)
        % ---------------------------------------------------------------
        
        % Change directory
        cd ../Source/

        % Grid perturbation coefficient
        PG = 0.99;

        % Generate nonuniform grid 
        [xx,yy,~,~,~,~,~,xx1,yy1,M] = GridGenerator(0,2,0,1,10,10,PG); 
        
        % Change directory
        cd(Directory)
    
        % Plot grid
        PlotGrid(xx,yy,M,xx1,yy1,'b','k')

    % ====================================================================
    %                            Figure 4
    % ====================================================================

    elseif strcmp(nfig,'4(a)')

        % ---------------------------------------------------------------
        %                           Figure 4(a)
        % ---------------------------------------------------------------

        % Plot L-shape domain grid
        PlotComplexShapes('L-shape',0)

    elseif strcmp(nfig,'4(b)')

        % ---------------------------------------------------------------
        %                           Figure 4(b)
        % ---------------------------------------------------------------

        % Plot circular domain grid
        PlotComplexShapes('Circle',0)

    elseif strcmp(nfig,'4(c)')

        % ---------------------------------------------------------------
        %                           Figure 4(c)
        % ---------------------------------------------------------------

        % Plot square domain with circular hole grid
        PlotComplexShapes('SquareWithHole',0)

    elseif strcmp(nfig,'4(d)')

        % ---------------------------------------------------------------
        %                           Figure 4(d)
        % ---------------------------------------------------------------

        % Plot L-shape domain grid with bonds
        PlotComplexShapes('L-shape',1)

    elseif strcmp(nfig,'4(e)')

        % ---------------------------------------------------------------
        %                           Figure 4(e)
        % ---------------------------------------------------------------

        % Plot circular domain grid with bonds
        PlotComplexShapes('Circle',1)

    elseif strcmp(nfig,'4(f)')

        % ---------------------------------------------------------------
        %                           Figure 4(f)
        % ---------------------------------------------------------------

        % Plot square domain with circular hole grid with bonds
        PlotComplexShapes('SquareWithHole',1)

    % ====================================================================
    %                            Figure 5
    % ====================================================================

    elseif strcmp(nfig,'5(a)')

        % ---------------------------------------------------------------
        %                           Figure 5(a)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for uniform grid & FA algorithm
        PlotNeighborList(1)

    elseif strcmp(nfig,'5(b)')

        % ---------------------------------------------------------------
        %                           Figure 5(b)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for uniform grid & PA-AC algorithm
        PlotNeighborList(2)

    elseif strcmp(nfig,'5(c)')

        % ---------------------------------------------------------------
        %                           Figure 5(c)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for uniform grid & IPA-AC algorithm
        PlotNeighborList(3)

    elseif strcmp(nfig,'5(d)')

        % ---------------------------------------------------------------
        %                           Figure 5(d)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for regular grid & FA algorithm
        PlotNeighborList(4)

    elseif strcmp(nfig,'5(e)')

        % ---------------------------------------------------------------
        %                           Figure 5(e)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for perturbed regular grid & FA algorithm
        PlotNeighborList(5)

    elseif strcmp(nfig,'5(f)')

        % ---------------------------------------------------------------
        %                           Figure 5(f)
        % ---------------------------------------------------------------

        % Plot discretized domain with all bonds for irregular grid & FA algorithm
        PlotNeighborList(6)

    % ====================================================================
    %                            Figure 6
    % ====================================================================

    elseif strcmp(nfig,'6(a)')

        % ---------------------------------------------------------------
        %                           Figure 6(a)
        % ---------------------------------------------------------------

        % Plot family and neighbor areas for the FA algorithm
        PlotNeighborAreasCoords(3.1,'FA')

    elseif strcmp(nfig,'6(b)')

        % ---------------------------------------------------------------
        %                           Figure 6(b)
        % ---------------------------------------------------------------

        % Plot family and neighbor areas for the PA-AC algorithm
        PlotNeighborAreasCoords(3.1,'PA-AC')

    elseif strcmp(nfig,'6(c)')

        % ---------------------------------------------------------------
        %                           Figure 6(c)
        % ---------------------------------------------------------------

        % Plot family, neighbor areas, and centroids for the IPA-AC algorithm
        PlotNeighborAreasCoords(3.1,'IPA-AC')

    % ====================================================================
    %                            Figure 7
    % ====================================================================
        
    elseif strcmp(nfig,'7') 

        % Plot influence functions
        PlotInfluenceFunction(1,1.5,0.01,0.1,[0 0.5 1 3 5 7])

    % ====================================================================
    %                            Figure 8
    % ====================================================================

    elseif strcmp(nfig,'8(a)') 

        % ---------------------------------------------------------------
        %                           Figure 8(a)
        % ---------------------------------------------------------------    

        % Plot horizontal pre-notch with intact bonds
        PlotPreNotch(1,0)

    elseif strcmp(nfig,'8(b)') 

        % ---------------------------------------------------------------
        %                           Figure 8(b)
        % ---------------------------------------------------------------    

        % Plot inclined pre-notch with intact bonds
        PlotPreNotch(2,0)

    elseif strcmp(nfig,'8(c)')

        % ---------------------------------------------------------------
        %                           Figure 8(c)
        % ---------------------------------------------------------------

        % Plot random pre-notches with intact bonds
        PlotPreNotch(3,0)

    elseif strcmp(nfig,'8(d)')

        % ---------------------------------------------------------------
        %                           Figure 8(d)
        % ---------------------------------------------------------------

        % Plot horizontal pre-notch with intact and broken bonds
        PlotPreNotch(1,1)

    elseif strcmp(nfig,'8(e)')

        % ---------------------------------------------------------------
        %                           Figure 8(e)
        % ---------------------------------------------------------------

        % Plot inclined pre-notch with intact and broken bonds
        PlotPreNotch(2,1)

    elseif strcmp(nfig,'8(f)')

        % ---------------------------------------------------------------
        %                           Figure 8(f)
        % ---------------------------------------------------------------

        % Plot random pre-notches with intact and broken bonds
        PlotPreNotch(3,1)

    % ====================================================================
    %                            Figure 9
    % ====================================================================

    elseif strcmp(nfig,'9(a)')

        % ---------------------------------------------------------------
        %                           Figure 9(a)
        % ---------------------------------------------------------------

        % Plot collinear segments
        PlotSegmentsIntersection(1)

    elseif strcmp(nfig,'9(b)')

        % ---------------------------------------------------------------
        %                           Figure 9(b)
        % ---------------------------------------------------------------

        % Plot parallel (not collinear) segments
        PlotSegmentsIntersection(2)

    elseif strcmp(nfig,'9(c)')

        % ---------------------------------------------------------------
        %                           Figure 9(c)
        % ---------------------------------------------------------------

        % Plot non-parallel segments
        PlotSegmentsIntersection(3)

    elseif strcmp(nfig,'9(d)')

        % ---------------------------------------------------------------
        %                           Figure 9(d)
        % ---------------------------------------------------------------

        % Plot non-parallel segments with endpoint intersection
        PlotSegmentsIntersection(4)

    % ====================================================================
    %                            Figure 10
    % ====================================================================

    elseif strcmp(nfig,'10(a)')

        % ---------------------------------------------------------------
        %                           Figure 10(a)
        % ---------------------------------------------------------------

        % Plot undeformed grid to illustrate bond breaking
        PlotBondBreaking(0,0)

    elseif strcmp(nfig,'10(b)')

        % ---------------------------------------------------------------
        %                           Figure 10(b)
        % ---------------------------------------------------------------

        % Plot deformed grid with broken bonds to illustrate bond breaking
        PlotBondBreaking(1,0)

    elseif strcmp(nfig,'10(c)')

        % ---------------------------------------------------------------
        %                           Figure 10(c)
        % ---------------------------------------------------------------

        % Plot deformed grid with fewer broken bonds due to a no-fail zone
        % to illustrate bond breaking
        PlotBondBreaking(1,1)

    % ====================================================================
    %                            Figure G4
    % ====================================================================

    elseif strcmp(nfig,'G4')

        % Plot discrete system for computation of fracture energy
        PlotNumericalFractureEnergy(1,6)

    else

        error('Figure identifier unknown.')

    end

end
