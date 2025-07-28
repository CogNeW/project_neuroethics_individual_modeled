function [params, negll, bic] = get_mtDDM_params(values, choices, rts)
% [params, negll, bic] = get_mtDDM_params(values, choices, rts)
%
% Estimate mtDDM parameters using genetic algorithms
%
% inputs:
% values: real value differences (top - bottom boundary) experienced by the
%     participant.
%     size nTrials x 2
%     values are positive if in favor of the option coded as choice=1
%     values are negative if in favor of the option coded as choice=-1
% choices: real choices a participant made. size nTrials x 1
%     1 or -1. 1 should be a choice in favor of the top boundary, -1 for
%    the bottom boundary.
% rts: real response times (in cs). size nTrials x 1
% 
% outputs:
% params: best estimate from grid search. Columns:
%   [drift slope for attribute 1, drift slope for attribute 2, 
%   drift latency 1, drift latency 2, barrier/boundary size]
% negll: negative log likelihood of the best-fitting grid search
% bic: bayesian information criterion value
%
% NOTE!!! input rts need to be in cs. output slopes and latencies are in cs - 
%     multiply them by 10 to bring into milliseconds
%
% MATLAB requirements are:
% Optimization Toolbox
% Global Optimization Toolbox
% MATLAB R2014b or higher
% 
% v1.2 Written by Nikki Sullivan and available as part of open access code for:
% Sullivan, NJ, Huettel, SA (2021) Healthy choice depends on the latency 
% and rate of information accumulation. Nature Human Behavior.
% v1.2 options added to GA to optimize search

    % constrain estimation for speed
    % NOTES: 
    % 1) These may need to change depending on input values.
    % 2) You can omit these by leaving them empty (e.g. lb = [];) but this 
    %   will take a long time to estimate.
    % 3) Note that the upper and lower estimation limits below for the
    %   boundary size is the boundary's absolute value
    lb = [0    0  1         1         .5];
    ub = [.02 .02 mean(rts) mean(rts) 2];

    options.Display = 'off';
    options.InitialPopulationRange = [lb;ub]; % where to start with parents
    % these last two slow things down but improve performance
    options.PopulationSize = 300; 
    options.MaxGenerations = size(lb,2)*200;
    
%     Testing options so that it runs quickly
    options.PopulationSize = 4;
    options.MaxGenerations = 10;

    [params, negll] = ga(@(x)(llFun(values, choices, rts, x)), ...
        size(lb,2),[],[],[],[], lb, ub, [], [3 4], options);

    bic = -2 .* (-negll) + size(lb,2) .* log(length(choices));

end

function negll = llFun(vals, dataChoices, dataRTs, params)
% negll = llFun(vals, dataRTs, dataChoices, params)
%
% calculates negative log likelihood for this set of parameters
%
% inputs:
%   vals: value differences, n x 2 where n is number of trials
%   dataRTs: RTs in cs, n x 1 where n is number of trials
%   dataChoices: choices (1 or -1), n x 1 where n is number of trials
%   params: ddm parameters, [slope1 slope2 lat1 lat2 bounds]
%           lat1 and lat2 are indices so must be integers
%
% note: vals and dataChoices should be calculated with 
%   respect to each other. For example if the value difference is
%   calcuated as left-right option, dataChoices should be 1 for the left
%   option and -1 for the right option.

    nTrials = length(dataChoices);

    % NOTE: below five RT bins are used. You may need to change this
    % depending on how disperse your RT distributions are. I recommend
    % doing some simulations and recoveries that match the idiosyncrasies
    % of your dataset to see what's best. Alternatively you can try not
    % specifying the number of bins and allow histcounts to choose this for
    % you; this did not work as well for the data tested here, however.
    % bin real RTs (log-spaced, as RT distributions are positively skewed)
    [~, RTBins, ~] = histcounts(log(dataRTs),5);
    RTBins = exp(RTBins); % log(cs) --> cs

    % how many simulations of each value difference pair?
    % 1000 is good, 5000 is better, and simulations indicate diminishing 
    % marginal returns after 10000 but 10000 is very slow.
    nSims = 5000; 
    
    nSteps = ceil(RTBins(end)); % max time steps allowed for simulations
    possDiffs = unique(vals, 'rows'); % get unique value difference pairs

    rvsNoise = .1; % standard deviation of noise distribution
    bound = [params(end) -params(end)]; % here's where you could add asymmetrical bounds
    dRates = [params(1)*possDiffs(:,1) params(2)*possDiffs(:,2)];

    % which latency is first? if equal, this is equivalent to an mDDM
    ndt = params([3 4]);
    [~, first] = min(ndt);
    [~, second] = max(ndt);

    % get value signal for each unique value difference pair
    % rows = value differences; columns = time steps
    rvs = [zeros(size(possDiffs,1), ndt(first)-1) ...
        repmat(dRates(:,first), 1, ndt(second)-ndt(first)) ...
        repmat(sum(dRates,2), 1, nSteps-ndt(second)+1)];
    % rvs(:,1) = ?; % FYI here's where you would add a starting point bias

    
    % simulate choices and RTs with this parameter set
    rt = NaN(nSims, 1);
    choice = rt;
    topChoice = false(nSims, 1);
    botChoice = topChoice;
    prob = zeros(2, length(RTBins), length(possDiffs));
    for d = 1:size(possDiffs,1)

        rvsPaths = repmat(rvs(d,:),nSims,1) + ... % paths
            normrnd(0, rvsNoise, nSims, nSteps);  % noise
        rvsSum = cumsum(rvsPaths,2); % create cumulative signal

        % reached top boundary
        [reachedT, timeT] = max(rvsSum >= bound(1), [], 2);
        % reached bottom boundary
        [reachedB, timeB] = max(rvsSum <= bound(2), [], 2);

        topChoice(:) = (reachedT & ~reachedB) | ... % got to top and never also bottom
            (reachedT & reachedB & (timeT < timeB)); % or reached top first
        botChoice(:) = (~reachedT & reachedB) | ... % got to bottom and never also top
            (reachedT & reachedB & (timeT > timeB)); % or reached bot first

        choice(topChoice) = 1;
        choice(botChoice) = -1;    
        rt(topChoice) =  timeT(topChoice);
        rt(botChoice) =  timeB(botChoice);    

        % which bin do the simulated RTs and choices fall into?
        [~, ~, whichRTBin] = histcounts(rt, RTBins);

        % number of responses in each response and RT bin
        for t = unique(whichRTBin(whichRTBin>0))' % for experienced RT bins
            prob(1:2, t, d) = ...
                [sum(choice==1 & whichRTBin==t) ... % top barrier choice
                sum(choice==-1 & whichRTBin==t)];   % bottom barrier choice
        end

        % simulations indicate doing this instead of clearing each iteration is faster
        rt(:) = NaN; 
        choice(:) = NaN;

    end
    % number -> proportion
    prob = prob ./ nSims;


    % evaluate fit to real data
    % find which bin real responses in:
    alltheRespBins = NaN(1,nTrials);
    alltheRespBins(dataChoices == 1) = 1;
    alltheRespBins(dataChoices == -1) = 2;
    % find which bin real RTs are in:
    [~, ~, allTheRTBins] = histcounts(dataRTs, RTBins);
    % for each real trial, find its likelihood given these params
    like = zeros(nTrials, 1);
    for trial = 1:nTrials

        % find real trial's value bin:
        valBin = possDiffs(:, 1) == vals(trial, 1) ...
            & possDiffs(:, 2) == vals(trial, 2);

        % assign the probability
        like(trial) = prob(alltheRespBins(trial), ...
            allTheRTBins(trial), valBin);

    end
    negll = -nansum(log(like+eps));

end