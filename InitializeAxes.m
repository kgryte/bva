function handles = InitializeAxes(handles)
% INITIALIZEAXES
%   Initialize axes within a figure window according to defined
%   specifications. Three subplots are generated: [i] 2D histogram subplot
%   and [ii] 2 1D histogram suplots along the upper and right edges of the
%   2D histogram subplot. Orientations of the 1D histogram subplots are
%   such that they are the projections of the data along the abscissa and
%   ordinate axes, respectively.
%
%   INPUTS:
%       HANDLES: structure with the following fields:
%          FIGURE: handle to the figure window. Suggested figure size:
%          [520, 196, 900, 700].
%           
%
%   OUTPUTS:
%       HANDLES: structure with the following fields:
%           SUBPLOT: nested structure with the following fields:
%               BVA: nested structure with the following fields:
%                   AXES: handles to the 2D Histogram axes (subplot).
%           SUBPLOT: nested structure with the following fields:
%               HISTOGRAM: nested structure with the following fields:
%                   MEAN: nested structure with the following fields:
%                      AXES: handle to the axes into which the histogram
%                         will be plotted.
%                   STDEV: nested structure with the following fields:
%                      AXES: handle to the axes into which the histogram
%                         will be plotted.
%
%
%
%   DEPENDENCIES:
%       (none)
%
%   -----------------------
%   EXAMPLE USAGE:
%       handles = InitializeAxes(handles);
%
%   ------------------------
%   History: 
%       [1] Kristofer Gryte. University of Oxford (2011).
%
%   Copyright (c) 2011. Kristofer Gryte. University of Oxford.
%
%   
%
%   Version: 1.0 (2011-08-02)
%
%   ----------------------
%   Author Information:
%       Kristofer Gryte
%       Clarendon Laboratory
%       University of Oxford
%       Oxford, United Kingdom
%       e-mail: k.gryte1 (at) physics.ox.ac.uk
%
%   ------------------------
%   TODO:
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

% Activate figure window:
figure(handles.figure);

% Configure the figure colormap:
load('BVAcolormap.mat', 'BVAcolormap');

set(handles.figure,...
    'Colormap', BVAcolormap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Subplots:

% Mean Histogram Subplot:
handles.Subplot.Histogram.Mean.Axes = subplot(3,1,1);

% Standard Deviation Histogram Subplot:
handles.Subplot.Histogram.StDev.Axes = subplot(3,1,2);

% BVA Subplot:
handles.Subplot.BVA.Axes = subplot(3,1,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manipulate Subplots:

% Get the current units:
OldUnits = get(handles.Subplot.Histogram.Mean.Axes, 'Units');

% Set Units: % Note: Change to characters...as the number of pixels is
% relative to user monitor...thus, standardize...
set([handles.Subplot.Histogram.Mean.Axes,...
    handles.Subplot.Histogram.StDev.Axes,...
    handles.Subplot.BVA.Axes],...
    'Units', 'pixels');

% Set Mean Histogram Position:
set(handles.Subplot.Histogram.Mean.Axes,...
    'Position', [150 500 350 150]);

% Set StDev Histogram Position:
set(handles.Subplot.Histogram.StDev.Axes,...
    'Position', [600 95 150 350]);

% Set BVA Subplot Position:
set(handles.Subplot.BVA.Axes,...
    'Position', [150 95 350 350]);

% Return units to 'characters':
set([handles.Subplot.Histogram.Mean.Axes,...
    handles.Subplot.Histogram.StDev.Axes,...
    handles.Subplot.BVA.Axes],...
    'Units', OldUnits); % Maybe Normalize?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF

