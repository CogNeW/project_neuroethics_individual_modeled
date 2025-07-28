function exampleScriptToEstimate_mtDDM
% exampleScriptToEstimate_mtDDM
%
% File to illustrate how to run the mtDDM estimation scripts

%% information on input data
%
% Choices
% choices need to be 1 or -1 and are n x 1 with n being the number of trials
% for the example data here, 1 is left choice, -1 is right choice
%
% Value inputs - these are the attribute values (e.g. taste and health)
% values are n x 2 with n being the number of trials
% values are positive if in favor of the option coded as choice=1
% values are negative if in favor of the option coded as choice=-1
% for the example data here, values are the left - right option value 
% differences for [taste health] 
%
% Response times
% rts are in centaseconds and are n x 1 with n being the number of trials

% load sample values, choices, and rts
load exampleData.mat


%% run estimation

% params variable is 1 x 5 in the following order: 
%   [drift slope for attribute 1, drift slope for attribute 2, 
%   drift latency 1, drift latency 2, barrier/boundary size]
%
%   units for slopes are units/cs, for latencies are cs.

[params, negll, bic] = get_mtDDM_params(values, choices, rts);

save('testDDM.mat','params','negll','bic')


end