function [DATA, REMOVEBURSTS] = FilterData(...
    ArrivalMatrix, FILTERS)
% FILTERDATA
%   Filter a photon arrival matrix using supplied filter parameters.
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
%       burst:
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
%       REMOVEBURSTS: a vector providing indices for the bursts which
%       failed filter criteria.
%
%   DEPENDENCIES:
%       SPLITVEC.M
%
%   -----------------------
%   EXAMPLE USAGE:
%       [DATA, REMOVEBURSTS] = FilterData(ArrivalMatrix, FILTERS);
%
%   ------------------------
%   History: 
%       [1] Kristofer Gryte. University of Oxford (2011).
%
%   Copyright (c) 2011. Kristofer Gryte. University of Oxford.
%
%   
%
%   Version: 1.0 (2011-07-30)
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
%       [1] Add checks for input variables.
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

REMOVEBURSTS = [];  

Apply = true; % FLAG to determine if apply filters or not. THIS IS NEEDED, FOR NOW, EVEN IF NO FILTERS ARE APPLIED!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Single Data Matrix:

% Replace all the NaNs:
ArrivalMatrix(isnan(ArrivalMatrix(:,3)),3) = 0;

% Return the assigned IDs of the various runs: % recall that non-burst
% photons are labeled with a '0'!!!!
FLAGVALS = SplitVec(...
    ArrivalMatrix(:,3),...ArrivalMatrix(~isnan(ArrivalMatrix(:,3)),3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'firstval'); % return just the first index marking a run.


% --- %
% [1] Start Times:
STARTINDICES = SplitVec(...
    ArrivalMatrix(:,3),...ArrivalMatrix(~isnan(ArrivalMatrix(:,3)),3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'first'); % return just the first index marking a run.

STARTINDICES(FLAGVALS == 0) = []; % NOTE: 

DATA.StartTimes = ArrivalMatrix(STARTINDICES, 1);

% --- %
% [2] End Times:
ENDINDICES = SplitVec(...
    ArrivalMatrix(:,3),...ArrivalMatrix(~isnan(ArrivalMatrix(:,3)),3),... % exclude any NaNs
    'equal',... % find all those runs of the same index (burst ID)
    'last'); % return just the last index completing a run.

ENDINDICES(FLAGVALS == 0) = []; % NOTE: 

DATA.EndTimes = ArrivalMatrix(ENDINDICES, 1);

% --- %
% [3] Duration:
DATA.Duration = DATA.EndTimes - DATA.StartTimes;

% --- %
% [4] Length:
DATA.Length = ENDINDICES - STARTINDICES;

% --- %
BURSTCELL(:,1) = SplitVec(...
    ArrivalMatrix(:, 3),...ArrivalMatrix(~isnan(ArrivalMatrix(:,3)), 2),...
    'equal',...
    'loc'); % this is a cellular array!

BURSTCELL(FLAGVALS == 0) = []; % NOTE:

% ---------- %
% Free up some memory:
% clear ArrivalMatrix
% ---------- %

% [5] DexDem: % DexDem = 0
DATA.DexDem(:,1) = cellfun(@(x) sum(ArrivalMatrix(x,2) == 0), BURSTCELL, 'UniformOutput', true);

% [6] DexAem: % DexAem = 2
DATA.DexAem(:,1) = cellfun(@(x) sum(ArrivalMatrix(x,2) == 2), BURSTCELL, 'UniformOutput', true);

% [7] AexAem: % AexAem = 1
DATA.AexAem(:,1) = cellfun(@(x) sum(ArrivalMatrix(x,2) == 1), BURSTCELL, 'UniformOutput', true);

% --- %
% [8] Sum Photons: DexDem + DexAem
DATA.SumPhotons = DATA.DexDem + DATA.DexAem;

% --- %
% [9] Total Photons: DexDem + DexAem + AexAem
DATA.TotalPhotons = DATA.SumPhotons + DATA.AexAem;

% --- %
% [10] Efficiency:
DATA.Efficiency = DATA.DexAem ./ (DATA.SumPhotons);

% --- %
% [11] Stoichiometry:
DATA.Stoichiometry = DATA.SumPhotons ./ (DATA.TotalPhotons);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Threshold Data:

if Apply == true
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Duration:
    if ~isempty(FILTERS.Duration.Min)

        IDs = find(DATA.Duration < FILTERS.Duration.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Duration(DATA.Duration < FILTERS.Duration.Min,:) = NaN;        
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Duration:
    if ~isempty(FILTERS.Duration.Max)

        IDs = find(DATA.Duration > FILTERS.Duration.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Duration(DATA.Duration > FILTERS.Duration.Max,:) = NaN;
           
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Length:
    if ~isempty(FILTERS.Length.Min)

        IDs = find(DATA.Length < FILTERS.Length.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Length(DATA.Length < FILTERS.Length.Min,:) = NaN;        
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Length:
    if ~isempty(FILTERS.Length.Max)

        IDs = find(DATA.Length > FILTERS.Length.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Length(DATA.Length > FILTERS.Length.Max,:) = NaN;        
        
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Min DexcDem:
    if ~isempty(FILTERS.DexDem.Min)

        IDs = find(DATA.DexDem < FILTERS.DexDem.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];
        
        DATA.DexDem(DATA.DexDem < FILTERS.DexDem.Min,:) = NaN;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Max DexcDem:
    if ~isempty(FILTERS.DexDem.Max)

        IDs = find(DATA.DexDem > FILTERS.DexDem.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.DexDem(DATA.DexDem > FILTERS.DexDem.Max,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Min DexcAem:
    if ~isempty(FILTERS.DexAem.Min)

        IDs = find(DATA.DexAem < FILTERS.DexAem.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.DexAem(DATA.DexAem < FILTERS.DexAem.Min,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max DexcAem:
    if ~isempty(FILTERS.DexAem.Max)

        IDs = find(DATA.DexAem > FILTERS.DexAem.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.DexAem(DATA.DexAem > FILTERS.DexAem.Max,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Min AexcAem:
    if ~isempty(FILTERS.AexAem.Min)

        IDs = find(DATA.AexAem < FILTERS.AexAem.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.AexAem(DATA.AexAem < FILTERS.AexAem.Min,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max AexcAem:
    if ~isempty(FILTERS.AexAem.Max)

        IDs = find(DATA.AexAem > FILTERS.AexAem.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.AexAem(DATA.AexAem > FILTERS.AexAem.Max,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Sum DD+DA:
    if ~isempty(FILTERS.SumPhotons.Min)

        IDs = find(DATA.SumPhotons < FILTERS.SumPhotons.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];
       
        DATA.SumPhotons(DATA.SumPhotons < FILTERS.SumPhotons.Min,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Sum DD+DA:
    if ~isempty(FILTERS.SumPhotons.Max)

        IDs = find(DATA.SumPhotons > FILTERS.SumPhotons.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.SumPhotons(DATA.SumPhotons > FILTERS.SumPhotons.Max,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Total Signal:
    if ~isempty(FILTERS.TotalPhotons.Min)

        IDs = find(DATA.TotalPhotons < FILTERS.TotalPhotons.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.TotalPhotons(DATA.TotalPhotons < FILTERS.TotalPhotons.Min,:) = NaN;
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Total Signal:
    if ~isempty(FILTERS.TotalPhotons.Max)

        IDs = find(DATA.TotalPhotons > FILTERS.TotalPhotons.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.TotalPhotons(DATA.TotalPhotons > FILTERS.TotalPhotons.Max,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Efficiency:
    if ~isempty(FILTERS.Efficiency.Min)

        IDs = find(DATA.Efficiency < FILTERS.Efficiency.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Efficiency(DATA.Efficiency < FILTERS.Efficiency.Min,:) = NaN;
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Efficiency:
    if ~isempty(FILTERS.Efficiency.Max)

        IDs = find(DATA.Efficiency > FILTERS.Efficiency.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];
        
        DATA.Efficiency(DATA.Efficiency > FILTERS.Efficiency.Max,:) = NaN;
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Min Stoichiometry:
    if ~isempty(FILTERS.Stoichiometry.Min)

        IDs = find(DATA.Stoichiometry < FILTERS.Stoichiometry.Min); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Stoichiometry(DATA.Stoichiometry < FILTERS.Stoichiometry.Min,:) = NaN;
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Stoichiometry:
    if ~isempty(FILTERS.Stoichiometry.Max)

        IDs = find(DATA.Stoichiometry > FILTERS.Stoichiometry.Max); % Column vector
        REMOVEBURSTS = [REMOVEBURSTS; IDs];

        DATA.Stoichiometry(DATA.Stoichiometry > FILTERS.Stoichiometry.Max,:) = NaN;
                
    end

end % end IF


% Retain only the unique indices:
REMOVEBURSTS = unique(REMOVEBURSTS);




% FOR NOW...
return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF
