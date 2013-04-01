function [DATA] = Confocal2ColorALExSolnMultDataBVA_Analysis(...
    varargin)
% CONFOCAL2COLORALEXSOLNMULTDATABVA_ANALYSIS 
%   Burst Variance Analysis (BVA) of confocal (2 color) solution 
%   alternating laser excitation (ALEx) data.
% 
%   INPUTS:
%       ARRIVALMATRIX: input data on which to perform Burst Variance
%       Analysis. The input matrix should be Nx3 or 3xN, with the following
%       as vector inputs: (1) Photon Arrival Times; (2) Detection Channel
%       (DexDem = 0; DexAem = 2; AexAem = 1; AexDem = 3); (3) Burst
%       Classification (those photons belonging to a burst should have a
%       corresponding index to a unique burst ID; background photons should
%       be indexed as NaN): ARRIVALMATRIX = ArrivalTimes | ID | BurstNumber
%
%       WINDOW: size of the sliding window over which FRET will be
%       calculated. Default: 5.
%
%       CLUSTERBINS: number of bins over which to cluster bursts
%       characterized by their mean FRET value. Mean standard deviation
%       will be calculated across bursts belonging to a particular bin and
%       compared against a test statistic to determine significance (and
%       here, the presence of dynamics). Default: 20.
%
%       CLUSTERMINWINDOWS: threshold for the minimum number of total windows
%       within a cluster from all bursts of the cluster to be analyzed 
%       (i.e., we need to make sure we have enough data to be worth our 
%       while). Default: 50.
%
%       CONFIDENCEINTERVAL: (1-alpha) to determine statistical significance
%       of cluster results. Default: 0.999.
%
%       FILTERS: nested structure providing minimum and maximum thresholds
%       for the following data streams:
%           DURATION: the temporal length of a detected burst [seconds]
%           LENGTH: the number of photons within a detected burst [photons]
%           DEXDEM: the number of photons detected in the donor channel
%               upon donor excitation.
%           DEXAEM: the number of photons detected in the acceptor channel
%               upon donor excitation. 
%           AEXAEM: the number of photons detected in the acceptor channel
%               upon acceptor excitation.
%           SUMPHOTONS: the number of photons detected in the donor &
%               acceptor channels upon donor excitation.
%           TOTALPHOTONS: the number of photons detected in the donor &
%               acceptor channels upon donor and acceptor excitation.
%               (excluding photons detected in donor channel upon acceptor
%               excitation: AexDem)
%           EFFICIENCY: the FRET efficiency.
%           STOICHIOMETRY: the relative fluorophore brightness.
%
%   OUTPUTS:
%       DATA: structure with the following data streams for each detected
%       burst: (NOTE: data streams are filtered according to input FILTERS)
%           DURATION: the temporal length of a detected burst [seconds]
%           LENGTH: the number of photons within a detected burst [photons]
%           DEXDEM: the number of photons detected in the donor channel
%               upon donor excitation.
%           DEXAEM: the number of photons detected in the acceptor channel
%               upon donor excitation. 
%           AEXAEM: the number of photons detected in the acceptor channel
%               upon acceptor excitation.
%           SUMPHOTONS: the number of photons detected in the donor &
%               acceptor channels upon donor excitation.
%           TOTALPHOTONS: the number of photons detected in the donor &
%               acceptor channels upon donor and acceptor excitation.
%               (excluding photons detected in donor channel upon acceptor
%               excitation: AexDem)
%           EFFICIENCY: the FRET efficiency.
%           STOICHIOMETRY: the relative fluorophore brightness.
%           RESULTS: an Nx2 array, where N is the number of bursts (after
%               filtering). The first column corresponds to the mean Efficiency
%               value for a burst and the second to the burst's Efficiency 
%               standard deviation. Array: E(efficiency) | sqrt(V(efficiency))
%           CLUSTERS: nested structure with the following fields:
%               CENTERS: vector containing the centers of the individual
%                  clusters.
%               STDEVDISTRMEAN: the Monte Carlo simulated standard
%                  deviation prediction.
%               STATISTICS: an Mx2 array, where M is the number of clusters
%                   ('slices' along the efficiency axis, into which bursts are
%                   grouped). The first column corresponds to the mean
%                   Efficiency value for the cluster, and the second to the
%                   standard deviation of the burst's within that cluster.
%               WINDOWMEANS: cellular array where each entry contains the
%                   the Efficiency values for all sliding windows for bursts
%                   belonging to a cluster.
%               DISTANCEFROMUPPERCI: a vector of the distance between the
%                   mean standard deviations of each cluster from the upper
%                   confidence bound. This expresses something about the extent
%                   of the 'dynamics'.
%           CONFIDENCEINTERVAL: nested structure with the following fields:
%               UPPER: the upper confidence bound generated from the Monte
%                  Carlo simulation and multiple hypothesis test.
%               LOWER: the lower confidence bound generated from the Monte
%                  Carlo simulation and the multiple hypothesis test.
%
%   DEPENDENCIES:
%       SPLITVEC.M
%       INTERSECT_SORTED.M (Lightspeed)
%       UNION_SORTED_ROWS.M (Lightspeed)
%       ISMEMBER_SORTED.M (Lightspeed)
%       SETDIFF_SORTED.M (Lightspeed)
%       FILTERDATA.M (private)
%       MONTECARLOPREDICTION.M (private; requires STATS toolbox)
%
%   -----------------------
%   EXAMPLE USAGE:
%       Confocal2ColorALExSolnMultDataBVA_Analysis(ArrivalMatrix);
%
%       Confocal2ColorALExSolnMultDataBVA_Analysis(ArrivalMatrix, 'window',
%       5);
%
%       Confocal2ColorALExSolnMultDataBVA_Analysis(ArrivalMatrix, 'window',
%       5, 'clusterbins', 20, 'clusterminwindows', 50);
%
%   ------------------------
%   History: 
%       [1] Joseph Torella, Yusdi Santoso. University of Oxford (2009).
%       [2] Johannes Hohlbein. University of Oxford (2010).
%       [3] Kristofer Gryte. University of Oxford (2011).
%
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
%   Copyright (c) 2011. Kristofer Gryte. University of Oxford.
%
%   
%
%   Version: 1.0 (2011-07-29)
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
%       [1] Provide appropriate checks on filter values (e.g., ensure
%       numeric, positive, etc.)
%       [2] Make number of Monte Carlo simulations an input variable?
%       [3] Address the two possibilities: (i) burst is split, thus sharing
%       a photon, possibly, in which case EndTime_n = StartTime_n+1 & (ii)
%       burst is single photon, thus EndTime_n = StartTime_n. For now, this
%       distinction has little effect, as this is rare and insignificant,
%       but for sake of completeness... 
%       [4] Explicitly provide the Bonferroni correction; here, our input
%       alpha for the CI assumes the correction has already been
%       considered. In reality, we should make the CI explicitly dependent
%       on CLUSTERBINS. Of course, requiring the user to consider the alpha
%       value beforehand should ideally force the user consider its effects
%       and adjust accordingly. Also, allow the parameter-value pair to be
%       input as 'alpha', instead of just the CI.
%
%   ------------------------
%   See also 
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOP



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StartUp/CleanUp:

StartUp();

C = onCleanup(@() CleanUp());


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

% ---- %
% Establish default values:
% ---- %
WINDOW = 5; % size of the sliding window (non-overlapping) over which FRET will be calculated; hence, for every 5 photons, we calculate FRET
CLUSTERBINS = 20; % number of bins over which we will cluster bursts characterized by their mean FRET value
CLUSTERMINWINDOWS = 50; % threshold for the minimum number of total windows within a cluster over all bursts for the cluster to be analyzed (i.e., we need to make sure we have enough data to be worth our while)

CONFIDENCEINTERVAL = 0.999; % (1-alpha) CI level to determine if observed variance is statistically significant
NUMMONTECARLOSAMPLES = 500; % Number of Monte Carlo simulations to run for CI bound generation

FILTERS = struct(...
    'DexDem',       struct('Min', [], 'Max', []),...
    'DexAem',       struct('Min', [], 'Max', []),...
    'AexAem',       struct('Min', [], 'Max', []),...
    'AexDem',       struct('Min', [], 'Max', []),...
    'Efficiency',   struct('Min', [], 'Max', []),...
    'Stoichiometry',struct('Min', [], 'Max', []),...
    'SumPhotons',   struct('Min', [], 'Max', []),...
    'TotalPhotons', struct('Min', [], 'Max', []),...
    'Duration',     struct('Min', [], 'Max', []),...
    'Length',       struct('Min', [], 'Max', []));

DATA = [];

% ---- %
% Parse input arguments:
% ---- %

ArrivalMatrix = varargin{1};

% Check!!!
if ~isnumeric(ArrivalMatrix)
    % error out!
    errordlg(['ERROR:invalid input argument for BVA function. Please input ',...
        'a numeric array for the initial argument.'],...
        'ERROR:invalid input argument');
    return;
end % end IF

% Perform check for data input size and ensure dimensionality is correct
% (column matrix)
if size(ArrivalMatrix,1) < size(ArrivalMatrix, 2) % More columns than rows
    % Take the tranpose:
    ArrivalMatrix = ArrivalMatrix';

elseif size(ArrivalMatrix,2) ~= 3 % too many or too few columns 
    
    % Error out!
    errordlg(['ERROR:invalid input argument for BVA function. Data input ',...
        'should be either (Nx3) or (3xN).'],...
        'ERROR:invalid input argument');
    return;
    
end % end IF/ELSEIF

% Check number of input arguments:
if (mod(nargin, 2) - 1) ~= 0 % odd number of input arguments
    % Throw a fit and error out!
    
    errordlg(['ERROR:invalid number of input arguments for BVA function. ',...
        'Please enter parameter-value pairs in addition to the data input.'],...
        'ERROR:invalid number of input arguments');
    return;
    
end % end IF

% Loop through and assign parameter-value pairs:
for NumArg = 2 : 2 : nargin
    
    switch lower(char(varargin{NumArg}))
        
        case {'win', 'window', 'winsize'}
            % Update the WINDOW size:
            temp = varargin{NumArg + 1};
            
            % Check!
            if ~isnumeric(temp)|| temp <= 1
                % Error out:
                answer = questdlg(...
                    ['NOTICE:An invalid value was provided for the WINDOW paramter. ',...
                    'Please supply an integer value greater than 1. ',...
                    'Would you like to continue with a default WINDOW = 5, ',...
                    'or return to the command-line and supply an appropriate input value?'],...
                    'WARNING:invalid input value',...
                    'Continue', 'Cancel',...
                    'Cancel');
                
                % Parse the answer:
                switch lower(answer)
                    case 'continue'
                        % Move along with default WINDOW value
                    case 'cancel'
                        % Abort!
                        return;
                end % end SWITCH
                
            else
                % Assign over the input WINDOW size:
                WINDOW = temp;
            end % end IF/ELSE
            
        case {'clusterbins', 'bins', 'numbins', 'binnumber', 'binnum'}
            
            % Update the number of cluster bins:
            temp = varargin{NumArg + 1};
            
            % Check!
            if ~isnumeric(temp)|| temp <= 1
                % Error out:
                answer = questdlg(...
                    ['An invalid value was provided for the CLUSTERBINS paramter. ',...
                    'Please supply an integer value greater than 1. ',...
                    'Would you like to continue with a default CLUSTERBINS = 20, ',...
                    'or return to the command-line and supply an appropriate input value?'],...
                    'WARNING:invalid input value',...
                    'Continue', 'Cancel',...
                    'Cancel');
                
                % Parse the answer:
                switch lower(answer)
                    case 'continue'
                        % Move along with default CLUSTERBINS value
                    case 'cancel'
                        % Abort!
                        return;
                end % end SWITCH
                
            else
                % Assign over the input CLUSTERBINS value:
                CLUSTERBINS = temp;
            end % end IF/ELSE
            
        case {'clusterminwindows', 'minwindows'}
            
            % Update the minimum number of windows within a cluster:
            temp = varargin{NumArg + 1};
            
            % Check!
            if ~isnumeric(temp)|| temp <= 1
                % Error out:
                answer = questdlg(...
                    ['An invalid value was provided for the CLUSTERMINWNDOWS paramter. ',...
                    'Please supply an integer value greater than 1. ',...
                    'Would you like to continue with a default CLUSTERMINWINDOWS = 50, ',...
                    'or return to the command-line and supply an appropriate input value?'],...
                    'WARNING:invalid input value',...
                    'Continue', 'Cancel',...
                    'Cancel');
                
                % Parse the answer:
                switch lower(answer)
                    case 'continue'
                        % Move along with default CLUSTERMINWINDOWS value
                    case 'cancel'
                        % Abort!
                        return;
                end % end SWITCH
                
            else
                % Assign over the input CLUSTERMINWINDOWS:
                CLUSTERMINWINDOWS = temp;
            end % end IF/ELSE
            
        case {'ci', 'confidenceinterval', 'conf', 'confidence', 'confinterval'}
            
            % Update the confidence interval:
            temp = varargin{NumArg + 1};
            
            % Check!
            if ~isnumeric(temp)|| (temp >= 1 || temp <= 0)
                % Error out:
                answer = questdlg(...
                    ['An invalid value was provided for the CONFIDENCEINTERVAL paramter. ',...
                    'Please supply a numeric value greater than 0 and less than 1. ',...
                    'Would you like to continue with a default CONFIDENCEINTERVAL = 0.999, ',...
                    'or return to the command-line and supply an appropriate input value?'],...
                    'WARNING:invalid input value',...
                    'Continue', 'Cancel',...
                    'Cancel');
                
                % Parse the answer:
                switch lower(answer)
                    case 'continue'
                        % Move along with default CONFIDENCEINTERVAL value
                    case 'cancel'
                        % Abort!
                        return;
                end % end SWITCH
                
            else
                % Assign over the input CONFIDENCEINTERVAL:
                CONFIDENCEINTERVAL = temp;
            end % end IF/ELSE
            
        case {'filters', 'filter', 'thresh', 'thresholds', 'threshold', 'filt'}
            
            % Update the FILTERS:
            temp = varargin{NumArg + 1};
            
            % Check!
            if ~isstruct(temp)
                % Error out:
                answer = questdlg(...
                    ['An invalid value was provided for the FILTERS paramter. ',...
                    'Please supply an appropriate structure containing relevant fields. ',...
                    'Would you like to continue with a default FILTERS (none), ',...
                    'or return to the command-line and supply an appropriate input?'],...
                    'WARNING:invalid input value',...
                    'Continue', 'Cancel',...
                    'Cancel');
                
                % Parse the answer:
                switch lower(answer)
                    case 'continue'
                        % Move along with default FILTERS value
                    case 'cancel'
                        % Abort!
                        return;
                end % end SWITCH
                
            else
                % Assign over the input FILTERS:
                FILTERS = temp;
            end % end IF/ELSE
            
        otherwise
            % Throw a notice that the input argument was not recognized:
            msgbox(sprintf(...
                ['NOTICE:Unrecognized input argument: %s. \n',...
                    'Parameter-value pair ignored']),...
                char(varargin{NumArg}));
            
    end % end SWITCH
    
end % end FOR


% ---- %
% Calculate remaining parameters:
PHOTONCUTOFF = floor(WINDOW/2); % Remove the first few photons of a burst (allows analyzed data to be centered within burst)

CLUSTERWIDTH = (1-0)/CLUSTERBINS; % How wide are our bins in which we are clustering bursts (UNITS: none --> FRET efficiency); (0,1) defines the range of our FRET values.

% ---- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering:

% Filter the data:
[DATA, REMOVEBURSTS] = FilterData(ArrivalMatrix, FILTERS); % this is a function!!!!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst Variance Analysis:

% Start the timer...
tic

% Launch the waitbar:
h.Waitbar = waitbar(...
    0,...
    'Commencing BVA...please be patient');

% ------------------------------------------------------------------- %
%% Winnow the Data:

% Replace all the NaNs:
ArrivalMatrix(isnan(ArrivalMatrix(:,3)),3) = 0;

% Return the assigned IDs of the various runs: % recall that non-burst
% photons are labeled with a '0'!!!!
FLAGVALS = SplitVec(...
    ArrivalMatrix(:,3),...ArrivalMatrix(~isnan(ArrivalMatrix(:,3)),3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'firstval'); % return just the first index marking a run.

% Extract some information about the bursts from the ArrivalMatrix:
STARTINDICES = SplitVec(...
    ArrivalMatrix(:,3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'first'); % return just the first index marking a run.

ENDINDICES = SplitVec(...
    ArrivalMatrix(:,3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'last'); % return just the last index completing a run.

STARTINDICES(FLAGVALS == 0) = []; % NOTE: 
ENDINDICES(FLAGVALS == 0) = [];

% Get the start times and end times for each burst:
STARTTIMES = ArrivalMatrix(STARTINDICES, 1);
ENDTIMES = ArrivalMatrix(ENDINDICES, 1);

% Grab the total number of bursts:
TOTALBURSTS = numel(STARTINDICES);

% Create an index array:
BURSTINDICES(:,1) = 1 : TOTALBURSTS; % column vector

% Check!!!! Before binning photons, need to ensure that no StartTimes
% and EndTimes are equal (i.e., no burst ends when another begins), as
% this creates problems when creating the edges (eliminates a non-burst
% bin).
AdjacentBursts = intersect_sorted(...
    STARTTIMES,...
    ENDTIMES);

if isempty(AdjacentBursts) == false
    % Adjacent bursts have been found!!!
    fprintf('%d bursts were found to be immediately adjacent.\n',...
        numel(AdjacentBursts));

    % Cycle through each Adjacent Burst:
    for ADJBURST = 1 : numel(AdjacentBursts)

        % Find the burst:
        AdjBurstIndex = find(...
            STARTTIMES == AdjacentBursts(ADJBURST),...
            1,...
            'first');

        % Increase the arrival time of any photons which
        % arrived at the previous start time by a fractional amount:
        ArrivalMatrix(ArrivalMatrix(:,1) == STARTTIMES(AdjBurstIndex, 1), 1) =...
            ArrivalMatrix(ArrivalMatrix(:,1) == STARTTIMES(AdjBurstIndex, 1), 1)...
            +...
            1e-10; % should be less than NI card time resolution...

        % Increase the start time of the next burst by a fractional amount:
        STARTTIMES(AdjBurstIndex, 1) =...
            STARTTIMES(AdjBurstIndex, 1)...
            +...
            1e-10;

    end % end FOR


end % end IF


% Create Edges:
Edges = union_sorted_rows(...
    STARTTIMES(:,1),...
    ENDTIMES(:,1) + 1e-10);

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
ArrivalMatrix(TF == true,:) = [];
BinIndices(TF == true) = [];    


% ------------------------------------------------------------------- %
%% Burst Removal:

% Remove the filtered bursts:
TF = ismember_sorted((BinIndices + 1) ./ 2, REMOVEBURSTS);

ArrivalMatrix(TF == true, :) = [];
BinIndices(TF == true, :) = [];


% Now, we need to rid ourselves of all photons originating from red
% excitation: DexDem = 0; DexAem = 2; AexAem = 1; AexDem = 3; => all
% photons labeled 1 or 3... (NOTE: this also removes acceptor-only
% bursts!)
LogicalIndexArray = (ArrivalMatrix(:,2) == 1 | ArrivalMatrix(:,2) == 3);
VALS2 = unique(BinIndices(LogicalIndexArray, :)); % This gives all burst indices which contain acceptor excitation photons

ArrivalMatrix(LogicalIndexArray,:) = []; % All acceptor excitation photons removed.
BinIndices(LogicalIndexArray,:) = []; % Note: this array will be important now, as all photons can be grouped by their unique identification with a Bin...
% NOTE: now VALS and unique(BinIndices) do not necessarily
% correspond!!!

% Determine which, if any, bursts have been removed from further
% analysis:
REMOVEBURSTS = [REMOVEBURSTS; (setdiff_sorted(VALS2, unique(BinIndices)) + 1) ./ 2];

% REMOVEBURSTS now contains all bursts which either did not meet the
% filtering criteria above or were acceptor-only...

% Frequently, bursts have values which are not numbers; find them:
NaNBURSTS = find(isnan(DATA.Efficiency(:,1)));

REMOVEBURSTS = [REMOVEBURSTS; NaNBURSTS];





% ------------------------------------------------------------------- %   
% Clustering: cluster the bursts according to mean Efficiency values:

% Define the cluster bounds: % NOTE:this creates an extra bin!!!
CLUSTEREDGES = -CLUSTERWIDTH/2:CLUSTERWIDTH:1+CLUSTERWIDTH/2;
CLUSTEREDGES = CLUSTEREDGES';

% Each burst is histogrammed into a cluster bin:
[Counts, CLUSTERBINIDX] = histc(...
    DATA.Efficiency(:,1),...
    CLUSTEREDGES);



% ------------------------------------------------------------------- %

% Split the Arrival Matrix into Bursts:
IndexCell = SplitVec(...
    BinIndices,...
    'equal',... % Group all elements which are equal
    'loc'); % Return the locations of these elements in a cellular array

% Each cellular element of IndexCell corresponds to a burst; within each
% cell element, the individual photons are assigned a location
% corresponding to their index in the (filtered) ArrivalMatrix.


% ------------------------------------------------------------------- %
%% BVA:

waitbar(0,...
    h.Waitbar,...
    'Performing sliding window analysis...');

% Initialize a BurstVarianceAnalysis array:
DATA.Results = nan(TOTALBURSTS, 2);

% Loop through the bursts and perform BVA:
COUNTER = 0;
CLUSTERSTATS = cell(CLUSTERBINS+1, 1);

TRACKER = 0;
for BURST = 1 : TOTALBURSTS
    
    % Check!!!
    if ismember(BURST, REMOVEBURSTS) % Determine if BURST is blacklisted...
        % Blacklisted burst found! Move along to next burst...
        continue;
    else
        % Update our counter for our non-blacklisted bursts...
        COUNTER = COUNTER + 1;
    end % end IF

    % Grab the current (non-blacklisted) BURST indices:
    INDICES = IndexCell{COUNTER};

    % Remove a 'cutoff' number of photons:
    if numel(INDICES) < PHOTONCUTOFF
        % Go to next iteration...
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    else
        INDICES(1:PHOTONCUTOFF) = [];
    end % end IF/ELSE       


    % Calculate total number of photons:
    TOTALPHOTONS = numel(INDICES);        

    % Check!!!
    if TOTALPHOTONS < WINDOW
        % Go to next iteration, as not enough photons in 'burst'...
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    end % end IF       



    % We must apply a window to the burst photons; this window is
    % non-overlapping, and thus, in most cases, the number of windows
    % within a burst is not an integer, but rather, a remainder of
    % photons cannot be used.  Let us discard these photons, so we can
    % reshape our photon array: 
    RemPhotons = rem(TOTALPHOTONS, WINDOW);

    % We have some leftover photons; remove them...
    INDICES(end-RemPhotons+1:end) = [];


    % Reshape the photons from the arrival matrix into an array of
    % non-overlapping windows: (NOTE: photons within each window are
    % placed along the columns, and the window number within the burst
    % is placed along the rows)
    PhotonArray = transpose(reshape(ArrivalMatrix(INDICES,2), WINDOW, [])); % No apriori knowledge of how many windows within a burst

    % Calculate the Efficiency value (FRET) within each window: (NOTE:
    % we exploit the fact that DexDem photons are labeled 0):
    EfficiencyVector = sum((PhotonArray == 2), 2) ./ WINDOW; % where WINDOW = sum of DexDem + DexAem and the numerator is the number of DexAem photons
        
    % Calculate the standard deviation of the efficiency values:
    EfficiencyStDev = std(EfficiencyVector);

    % Place the mean and stdev of the efficiency in our BVA array:
    DATA.Results(BURST, :) = [mean(EfficiencyVector), EfficiencyStDev];  
    
    
    % Place the efficiency vector with the burst's respective cluster:
    WhichCluster = CLUSTERBINIDX(BURST);

    % Check!!!
    if WhichCluster == (CLUSTERBINS+2)
        % This accounts for values which are precisely equal to the
        % last edge.  We place this data in the last bin.
        WhichCluster = WhichCluster - 1;
    end % end IF/ELSE

    CLUSTERSTATS{WhichCluster} = [CLUSTERSTATS{WhichCluster}; EfficiencyVector];
    
    
    % Update our waitbar:
    waitbar(BURST/TOTALBURSTS, h.Waitbar);
    

end % end FOR    

fprintf(['%d bursts did not have sufficient photons ',...
    'to perform a sliding window calculation.\n\n'],...
    TRACKER);

% Tidy-up our Cluster Statistics:
RemoveClusters = cellfun(@isempty, CLUSTERSTATS);
CLUSTERSTATS(RemoveClusters) = [];
CLUSTEREDGES(RemoveClusters) = [];

RemoveClusters = cell2mat(cellfun(@(x) (numel(x)<CLUSTERMINWINDOWS),...
    CLUSTERSTATS,...
    'UniformOutput', false));
CLUSTERSTATS(RemoveClusters) = [];
CLUSTEREDGES(RemoveClusters) = [];


IDs = setdiff(BURSTINDICES, REMOVEBURSTS); % Retain only 'good' bursts

DATA.Results = DATA.Results(IDs,:);

% Remove any data in which the mean efficiency is either 0 or 1: (why?
% these 'events' tend to arise when having low photon statistics and will
% bunch up at the corner of the plots. Because of discrete photon
% statistics, these will occur often and affect the relatively weighting of
% the data elsewhere which has higher photon statistics and is more
% interesting (usually). Hence, we throw away this data, if present, and
% only keep that data which has more 'realistic' mean values.)
DATA.Results(DATA.Results(:,1) == 0 | DATA.Results(:,1) == 1, :) = [];

% Remove any data in which the standard deviation is zero: (why? Having a
% single STD is not realistic, as fluctuations should be present simply to
% to photon statistics. This also might be indicative of insufficient
% number of WINDOWS allowing for STD calculation. Once again, with low
% photon statistics, this is possible, and can bias the weighting of your
% 2D histogram.)
DATA.Results(DATA.Results(:,2) == 0, :) = [];

% Remove all NaN rows:
DATA.Results(isnan(DATA.Results(:,1)), :) = [];


% Tidy-up:
clear ArrivalMatrix
clear BinIndices

% ------------------------------------------------------------------- %
%% Cluster Statistics:

waitbar(0,...
    h.Waitbar,...
    'Analyzing cluster statistics...');


% ------------ %
% [ ] Group bursts with similar mean values and calculate population
% standard deviation:
ClusterCenters = CLUSTEREDGES(1:end-1) + CLUSTERWIDTH/2;

% Initialize some variables:
NumClusters = numel(CLUSTERSTATS);
HighCIvec = zeros(NumClusters, 1);
LowCIvec = zeros(NumClusters, 1);
StDevDistrMean = zeros(NumClusters, 1);

ClusterStats = nan(NumClusters, 2);

% Loop through each 'slice' (or 'cluster' of data) from the 2D Histogram
% and calculate the statistics:
for CLUSTER = 1 : NumClusters      

    ClusterStats(CLUSTER,2) = std(CLUSTERSTATS{CLUSTER});
    ClusterStats(CLUSTER,1) = mean(CLUSTERSTATS{CLUSTER});

    % Calculate StDev Distribution:
    Edges = 0 : (1-0)/1000 : 0.6; % why 0.6? Vertical upper bound (could be something else, but for practical purposes, we make this realistic)
    NumWindows = numel(CLUSTERSTATS{CLUSTER});
    StDevDistr = MonteCarloPrediction(...
        Edges,...
        NumWindows,... % number of windows
        WINDOW,... % number of photons within a window
        ClusterCenters(CLUSTER),... % probability of energy transfer
        NUMMONTECARLOSAMPLES); % number of samples generated by Monte Carlo simulation
    
    CumDistr = cumsum(StDevDistr);
    INV = 1 - CumDistr;
    StDevDistrMean(CLUSTER) = sum(StDevDistr .* Edges); % What is this? This should correspond to the theoretical shot noise prediction curve based on the number of photons specified by WINDOW

    %-- CI --
    %-- ADJUST CI BY THE  MULT HYPOTH TEST --
    absvec = abs(CumDistr - CONFIDENCEINTERVAL);
    [M, closest] = min(absvec);

    absvecLOW = abs(INV - CONFIDENCEINTERVAL);    
    closestLOW = find(absvecLOW == min(absvecLOW),...
        1, 'last');

    % Why are we calculating the CI off an index? The index tells us
    % something of about the number of false positives.
    HighCIvec(CLUSTER,:) = (closest-1)*.001; % 0.001 is 99.9% Confidence interval
    LowCIvec(CLUSTER,:) = (closestLOW-1)*.001;

    % Update our waitbar:
    waitbar(CLUSTER / NumClusters,...
        h.Waitbar);


end % end FOR

% Assign over to application data:
DATA.Clusters.Centers = ClusterCenters;
DATA.ConfidenceInterval.Upper = HighCIvec;
DATA.ConfidenceInterval.Lower = LowCIvec;
DATA.Clusters.StDevDistrMean = StDevDistrMean;
DATA.Clusters.Statistics = ClusterStats;
DATA.Clusters.WindowMeans = CLUSTERSTATS;
DATA.Clusters.DistanceFromUpperCI =...
    ClusterStats(:,2) - HighCIvec;


% ------------------------------------------------------------------- %
%% Tidy-up:


    
% Delete waitbar:
delete(h.Waitbar);


TOTALTIME = toc






return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StartUp:
function StartUp()
%
%
%

addpath(...
    'Other/lightspeed',...
    'Other/SplitVec');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CleanUp:
function CleanUp()
%
%
%

rmpath(...
    'Other/lightspeed',...
    'Other/SplitVec');



