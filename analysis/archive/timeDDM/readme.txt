-------------------------------------------------------------------------------
"Fast" mtDDM estimation files
-------------------------------------------------------------------------------
These estimation files can be used to run the mtDDM on your own datasets. For more on the mtDDM, see Sullivan, N.J., Huettel, S.A. (2021) 
Healthy choice depends on the latency and rate of information accumulation. Nature Human Behavior.


Files included:

exampleScriptToEstimate_mtDDM.m - this shows an example of loading a dataset and running the mtDDM estimation
get_mtDDM_params.m - this estimates mtDDM parameters
exampleData.mat - sample dataset.
get_mtDDM_params.R - an example of how to do the same thing using R instead of MATLAB.


-------------------------------------------------------------------------------
MATLAB requirements
-------------------------------------------------------------------------------
To run this code you will need:
MATLAB R2014b or higher
Optimization Toolbox
Global Optimization Toolbox
Deep Learning Toolbox


-------------------------------------------------------------------------------
Guidance
-------------------------------------------------------------------------------
the get_mtDDM_params function takes attribute values, choices, and response times from a participant and outputs a set of mtDDM parameters and estimates of how well those parameters fit participants' actual choices and response times: negative log likelihood values, and BIC (Bayesian Information Criterion). see "exampleScriptToEstimate_mtDDM.m" for an example of how to run this script.

The script is run with the following inputs and outputs:
[params, negll, bic] = get_mtDDM_params(values, choices, rts)

inputs
-------------------------------------------------------------------------------
It is important to ensure that choices and value difference pairs are coded with respect to each other. Details on this are below.

"values" variable: values needs to be n x 2 with n being the number of trials
This allows for two attributes (attribute 1, first column; attribute 2, second column) to guide choices and response times.
Values are positive if in favor of the option coded as choice=1
Values are negative if in favor of the option coded as choice=-1
For example, if choice==1 means the left choice, and choice==-1 means the right choice, then the value differences here should be encoded as "left minus right" such that larger positive numbers mean that the left option is more favored, and negative numbers mean the right option is more favored.
N.B., one source of very long estimation times is if the value differences are continuous rather than discreet or if the total number of unique value difference pairs is large. This is because the estimation script simulates choices and response times for a candidate set of parameters only for unique value pair combinations, which in a discrete case (e.g food ratings) is typically less than the total number of trials. It is recommended therefore to reduce the number of unique value pairs by, e.g., rounding or normalizing, if your input values are continuous (e.g. subjective values).
It also may make sense to normalize values such that both are within the same range (e.g. normalize both attribute values from -1 to 1) so that the magnitude of each drift slope can be compared.
If one value column is all zeros, this resolves to a simple DDM, though note that the BIC value should be modified as it is calculated assuming that there are five parameters, but if one attribute column is all zeros then there are essentially only three variables (one slope, one latency, one barrier).

"choices" variable: choices need to be 1 or -1 and are n x 1 with n being the number of trials

"rts" variable: response times for each choice, n x 1 with n being the number of trials
By convention, both input and output RTs are in centaseconds

!! Important:
This code constrains possible outcome parameters (see lines 44-45) which speeds up estimation. However, these upper and lower parameter constrains may well need to be different for your dataset depending, for example, on the range of your values. After estimation, it is best to check to ensure that participants are not fitting near any upper or lower constraint. Note that the boundary constraint is a constraint on a symmetrical boundary's absolute value.


Outputs:
-------------------------------------------------------------------------------

params output has the following columns: 
[drift slope for attribute 1, drift slope for attribute 2, drift latency 1, drift latency 2, barrier/boundary size]

Drift slopes set the amount of evidence accumulated for each attribute at each time step.
Drift latencies set the time at which each attribute begins being added to the value signal.
Barrier/boundary size sets the amount of evidence required before a choice is made.

The equation goes as follows:
V_t = V_t-1 + a_1 * d_1 * x_1 + a_2 * d_2 * x_2 + e_t, where 
a_1 = 0 before drift latency 1 is reached, and is 1 afterward;
a_2 = 0 before drift latency 2 is reached, and is 1 afterward

where V_t is the value signal at time point t, a is an attribute's latency indicator, d is an attribute's drift slope, x_1 is the value difference for attribute 1, x_2 is the same for attribute 2, and e is random i.i.d. Gaussian noise.


A note on interpretation of two drift slopes:
-------------------------------------------------------------------------------
It is theoretically reasonable (and mathematically equivalent) to formulate the mtDDM instead as:
V_t = V_t-1 + d * (a_1 * x_1 + a_2 * q * x_2) + e_t
where d is a temperature parameter and q is a discount on the value of attribute 2.

This formulation is convincingly argued for in Krajbich 2021 (https://psyarxiv.com/vnxsu/).

Converting the slope parameters here to a temperature and discount parameter is a simple matter of arithmetic. In fact, it may be more efficient to estimate these as independent slopes and then convert them to temperature and discount, as it can be more difficult to estimate the discount parameter if it has a very large range of possible values (e.g. from 0 to 10,000) and two slopes will have a smaller range in possible values. Determining which attribute receives the discount is, I suppose, up to you based on your theory as it doesn't matter either way mathematically.