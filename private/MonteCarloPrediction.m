function STDDISTRIBUTION = MonteCarloPrediction(EDGES, EVENTS,...
    TRIALS, PROBABILITY, ITERATIONS)
% MONTECARLOPREDICTION
%   Use Monte Carlo methods to generate a distribution of standard
%   deviation (std) values observed from a binomial process (Bernoulli
%   Trial) with a specified probability of success.
%
%   INPUTS:
%       EDGES: a vector defining the bin edges into which the simulated
%       standard deviations will be histogrammed.
%
%       EVENTS: the number of events (collections of trials) to be 
%       sampled from the distribution.
%
%       TRIALS: the number of trials.
%
%       PROBABILITY: the probability of success. For example, a series of
%       coin flips is a collection of events characterized by Bernoulli 
%       trials with probability 0.5.
%
%       ITERATIONS: the number of times to repeat each collection of events
%       so as to arrive at an average value for the calculated standard
%       deviations. Default: 500.
%
%
%   OUTPUTS:
%       STDDISTRIBUTION: a probability distribution of the standard
%       deviation with support defined by EDGES. Output is a vector.
%
%   DEPENDENCIES:
%       BINORND.M (STATS toolbox)
%
%   --------------------------
%   EXAMPLE USAGE:
%       MonteCarloPrediction([0:0.0001:1], 5, 100, 0.5);
%
%       MonteCarloPrediction([0:0.0001:1], 5, 100, 0.5, 500);
%
%   --------------------------
%   History: 
%       [1] Joseph Torella, Yusdi Santoso. University of Oxford (2009).
%       [2] Kristofer Gryte. University of Oxford (2011).
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
%   ---------------------------
%   TODO:
%       [1] Generate a function independent of the STATS toolbox.
%
%
%   ---------------------------
%   See also BINORND
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOP





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input Arguments:

% Check!
if nargin <= 4
    ITERATIONS = 500; % default value
end % end IF

% Check EDGES:
if size(EDGES,1) ~= 1 && size(EDGES,2) ~= 1
    % Throw an error:
    errordlg(['ERROR:invalid input argument. Please input a vector for the ',...
        'EDGES argument.'],...
        'ERROR:invalid input argument');
    
    return;
end % end IF

% Check EVENTS, TRIALS, PROBABILITY, and ITERATIONS
temp = [EVENTS, TRIALS, ITERATIONS];
if sum((~isinteger(temp)) + (temp>0)) ~= 6
    % Throw an error:
    errordlg(['ERROR:invalid input argument. Please input positive integers for the ',...
        'EVENTS, TRIALS, and ITERATIONS arguments.'],...
        'ERROR:invalid input argument');
    
    return;
end % end IF

if (PROBABILITY <= 0 || PROBABILITY > 1)
    % Throw an error:
    errordlg(['ERROR:invalid input argument. Please input a probability between 0 and 1 ',...
        'for the PROBABILITY argument.'],...
        'ERROR:invalid input argument');
    
    return;
end % end IF




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

% --- % 
% Parameters:
Runs = 3; % used in the event of a large number of requested standard deviations, in which case, for memory purposes, we need to break the work into smaller chunks.    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Monte Carlo Prediction:

if EVENTS*ITERATIONS <= 10000 % 10,000 ? A magical number discovered heuristically.
    
    % Generate the matrix of numbers drawn from a binomial distribution
    % with total trials equal to TRIALS and probability equal to
    % PROBABILITY:
    MyMatrix = binornd(...
        TRIALS,...
        PROBABILITY,...
        EVENTS,...
        ITERATIONS); % this is a matrix of size: EVENTS x ITERATIONS
    
    % Calculate the standard deviation along each column
    StdMatrix = std(MyMatrix ./ TRIALS);

elseif EVENTS*ITERATIONS > 10000
    
    % With many iterations, we run into the problem of very large matrices
    % and memory issues. Let us break it down more manageably:
    % ---- %
    
    RemainingSamples = ITERATIONS;
    
    % Pre-allocate our array:
    StdMatrix = nan(1, ITERATIONS); 
    
    % Initialize an index counter:
    Counter = 1;
    
    for k = 1 : Runs
        
        % Determine the number of samples to run:
        Samples = floor(RemainingSamples / (Runs-k+1));
        
        % Generate the matrix of numbers drawn from a binomial distribution
        % with total trials equal to TRIALS and probability equal to
        % PROBABILITY:
        MyMatrix = binornd(...
            TRIALS,...
            PROBABILITY,...
            EVENTS,...
            Samples);
        
        % Calculate the standard deviation along each column and append to
        % our standard deviation matrix:
        StdMatrix(Counter : (Counter+Samples-1)) = std(MyMatrix ./ TRIALS);
        
        % Calculate the remaining number of samples:
        RemainingSamples = RemainingSamples - Samples;
        
        % Update our counter:
        Counter = Counter + Samples;
            
    
    end % end FOR
    
       
end % end IF/ELSEIF
    
% Histogram the standard deviations where the bin edges are defined by
% EDGES:
[Counts, BinIndices] = histc(StdMatrix, EDGES);

% Determine the probability of a particular standard deviation, by
% normalizing by the total number of counts:
STDDISTRIBUTION = Counts / sum(Counts); % this is a vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF
