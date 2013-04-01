%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst Variance Analysis:
%
% Script:
%   This MATLAB script allows the user to perform Burst Variance Analysis
%   (BVA) and output this information to the user workspace and a figure.
%   Details concerning the BVA method and its application may be found in
%   the following publication:
%
%       Torella, JP, Holden, SJ, Santoso, Y, Hohlbein, J, Kapanidis, AN
%       (2011). Identifying molecular dyanmics in single-molecule FRET 
%       experiments with burst variance analysis. Biophys J.
%       100(6):1568-77.
%
%   To run the script, set the Current Directory to the folder containing
%   the BVA package and run the script BURSTVARIANCEANALYSIS.M either from
%   the editor or from the command-line using the following command:
%   ' BurstVarianceAnalysis '. 
%
%   The test data is real experimental data. The input data to be analyzed
%   must meet two requirements:
%
%       [1] Individual photons are sorted by their arrival times and are
%       sorted according to a ALEx (see Kapanidis) color identification 
%       scheme: (i) upon donor excitation (Dex): donor emission (Dem) -->
%       0, acceptor emission (Aem) --> 2; (ii) upon acceptor excitation
%       (Aex): donor emission (Dem) --> 3, acceptor emission (Aem) --> 1.
%       Thus, [DexDem, DexAem, AexAem, AexDem] --> [0, 2, 1, 3].
%
%       [2] A BURST SEARCH must have already been performed. Consult early
%       work out of the Siedel laboratory for this methodology. The two
%       most important outputs of the burst search are the start and end
%       times of the bursts. Such information can then be used to label
%       each photon as either belonging to a burst or not, and, if
%       belonging, to which burst it belongs.
%
%   Accordingly, the input data should be of the following format:
%
%       DATA = [ Arrival Times | Photon ID (0,2,1,3) | Burst ID]
%
%   As an example:
%       DATA = [ 0.40, 0, NaN;...
%                0.51, 1, NaN;...
%                0.52, 1,   1;...
%                0.53, 0,   1;...
%                0.54, 2,   1;...
%                0.55, 0,   1;...
%                0.62, 3, NaN;...
%                   .  .    .
%                   .  .    .
%                   .  .    .
%               51.70, 2, NaN;...
%               51.74, 0, 100;...
%               51.75, 2, 100;...
%               51.76, 2, 100;...
%               51.77, 1, 100;...
%               51.78, 1, 100;...
%               52.81, 2, NaN];
%
%   In the 'Testing' section, an example demonstrates how an Arrival Matrix
%   containing only arrival times and photon IDs can be transformed to
%   accommodate burst search input, as well as highlights some
%   considerations regarding the data (e.g., long bursts being split into
%   two, resulting in two adjacent bursts sharing the same photon, and the
%   need to increasing end times by a fractional amount to ensure
%   appropriate binning [NOTE: this, of course, can be avoided using
%   a binning method other than HISTC.M]).
%
%   Following the 'Testing' section, find the 'Initial Parameters' section,
%   which provides a starting point for input parameters. The BVA parameters
%   are typical values, as detailed in Torella, et al (2011); whereas, the
%   filter parameters will vary depending on the experimental setup and
%   conditions. For the example data included, the filter values are
%   found appropriate.
%
%   And finally, in the 'BVA' section, the figure generation, analysis, and
%   plotting is broken down into 4 steps (functions), with example usage
%   included.
%
%   For additional questions and concerns regarding the individual
%   functions themselves, please consult each function's HELP. Importantly,
%   more detail concerning the nature and structure of the input and output
%   parameters is provided. Further information concerning BVA
%   interpretation may be found in the documentation provided: HELP.pdf .
%   
%   Lastly, if any bugs or breaks are found in the code, please contact the
%   author of this package detailing the problem and provide any additional
%   evidence of the error/problem, which includes, but is not limited to,
%   error messages, figures, improper output data structure, et cetera. The
%   author's contact information may be found below.
%
%   NOTES:
%       [1] Need STATS toolbox for generating random numbers drawn from a
%       binomial distribution.
%
%   ----------------------
%   Publications: 
%       [1] Torella, JP, Holden, SJ, Santoso, Y, Hohlbein, J, Kapanidis, AN
%       (2011). Identifying molecular dyanmics in single-molecule FRET 
%       experiments with burst variance analysis. Biophys J. 100(6):1568-77.
%
%   NOTICE: 
%       If this software is used in the analysis of data which is to be
%       presented in publication, please cite the above work and
%       acknowledge the authors of this software.
%
%       For purposes of redistribution, use, and/or promotion of products
%       derived, consult the included license agreement.
%
%       Additionally, further license agreements are contained in
%       specialized subfunctions included in the directory 'Other'. Use,
%       redistribution, and derivative works are subject to any additional
%       requirements imposed therein.
%
%   Copyright (c) 2011. Kristofer Gryte. University of Oxford.
%
%   
%
%   Version: 1.0 (2011-08-04)
%
%   ----------------------
%   Author Information:
%       Kristofer Gryte
%       Clarendon Laboratory
%       University of Oxford
%       Oxford, United Kingdom
%       e-mail: k.gryte1 (at) physics.ox.ac.uk
%
%   ----------------------
%   TODO:
%       
%       [2] Documentation
%       
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing:

addpath(...
    'Other/lightspeed',...
    'Other/SplitVec');

clear ArrivalMatrix StartTimes EndTimes

load('MyArrivalMatrix.mat');

% Check!!!! Before binning photons, need to ensure that no StartTimes
% and EndTimes are equal (i.e., no burst ends when another begins), as
% this creates problems when creating the edges (eliminates a non-burst
% bin).
AdjacentBursts = intersect_sorted(...
    StartTimes,...
    EndTimes);

if isempty(AdjacentBursts) == false
    % Adjacent bursts have been found!!!
    fprintf('%d bursts were found to be immediately adjacent.\n',...
        numel(AdjacentBursts));

    % Cycle through each Adjacent Burst:
    for ADJBURST = 1 : numel(AdjacentBursts)

        % Find the burst:
        AdjBurstIndex = find(...
            EndTimes == AdjacentBursts(ADJBURST),...
            1,...
            'first');

        % Increase the arrival time of any photons which
        % arrived at the previous start time by a fractional amount:
        ArrivalMatrix(ArrivalMatrix(:,1) == EndTimes(AdjBurstIndex, 1), 1) =...
            ArrivalMatrix(ArrivalMatrix(:,1) == EndTimes(AdjBurstIndex, 1), 1)...
            +...
            1e-10; % should be less than NI card time resolution...

        % Increase the start time of the next burst by a fractional amount:
        EndTimes(AdjBurstIndex, 1) =...
            EndTimes(AdjBurstIndex, 1)...
            +...
            1e-10;

    end % end FOR


end % end IF


% Configure the ArrivalMatrix:
Edges = union_sorted_rows(...
    StartTimes,...
    EndTimes + 1e-10); % Bump the end times by a fractional amount so that last photon can be placed in appropriate burst bin.

% Histogram the arrival times:
[Counts, BinIndices] = histc(ArrivalMatrix(:,1), Edges);

% Determine which bins to throw away: (e.g., all those bins containing
% photons which do not belong to a burst)

% Discard the even BinIndices, as these correspond to bins defined by
% EdgeA = EndTime(1); EdgeB = StartTime(2). Hence, this bin (EdgeA,
% EdgeB) corresponds to photons which do not belong to a burst (i.e., photons
% occurring after the end of the previous burst and before the start of
% the next burst):
VALS1 = 0:2:max(BinIndices); % Note: we include 0 to account for photons not falling into any defined bin.
TF = ismember_sorted(BinIndices, VALS1);


% Remove all photons not belonging to a burst, thus creating a 'burst
% arrival matrix':
BinIndices(TF == true) = NaN; %[];

% Tack on the BinIndices (i.e., Burst Indices) to the Arrival Matrix:
ArrivalMatrix(:,3) = BinIndices;

clear BinIndices TF VALS1 Counts Edges


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Parameters:

WINDOW = 5;
CLUSTERBINS = 20;
CLUSTERMINWINDOWS = 50;
CONFIDENCEINTERVAL = 0.999;

% ---------- %
%  Filters:
% ---------- %

% DexDem:
FILTERS.DexDem.Min = [];
FILTERS.DexDem.Max = [];

% DexAem:
FILTERS.DexAem.Min = [];
FILTERS.DexAem.Max = [];

% AexAem:
FILTERS.AexAem.Min = [50];
FILTERS.AexAem.Max = [];

% Sum Photons: DexDem + DexAem
FILTERS.SumPhotons.Min = [50];
FILTERS.SumPhotons.Max = [];

% Total Photons: DexDem + DexAem + AexAem
FILTERS.TotalPhotons.Min = [];
FILTERS.TotalPhotons.Max = [];

% Efficiency:
FILTERS.Efficiency.Min = [];
FILTERS.Efficiency.Max = [];

% Stoichiometry:
FILTERS.Stoichiometry.Min = [0.45];
FILTERS.Stoichiometry.Max = [];

% Duration: [seconds]
FILTERS.Duration.Min = [];
FILTERS.Duration.Max = [];

% Length: [photons]
FILTERS.Length.Min = [100];
FILTERS.Length.Max = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BVA:

% [1] Create a figure:
handles.figure = figure;

set(handles.figure,...
    'Position', [520, 196, 900, 700]); % [x-pos, y-pos, width, height] % NOTE: the specified values are appropriate for the INITIALIZEAXES function to follow. Should the figure window be too large or small, subsequent re-adjustments in the INITIALIZEAXES function may be needed.

% [2] Initialize the figure and axes:
handles = InitializeAxes(handles);

% [3] Perform BVA:
[DATA] = Confocal2ColorALExSolnMultDataBVA_Analysis(...
    ArrivalMatrix,...
    'window', WINDOW,...
    'clusterbins', CLUSTERBINS,...
    'clusterminwindows', CLUSTERMINWINDOWS,...
    'confidenceinterval', CONFIDENCEINTERVAL,...
    'filters', FILTERS);

% [4] Plot the results:
handles = PlotResults(DATA, WINDOW, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOS






