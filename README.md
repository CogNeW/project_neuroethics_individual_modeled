# Code for 'The Dynamics of Neuroethical Decisions Depend on Psychological Traits'

> **Authors:**  
> Brian Kim<sup>1</sup>, Kiante Fernandez<sup>1,2</sup>, John D. Medaglia<sup>1,3</sup>  
> <sup>1</sup>Department of Applied Cognitive and Brain Sciences, Drexel University, Philadelphia, PA, 19104, USA  
> <sup>2</sup>Department of Psychology, University of California, Los Angeles, CA 90095, USA  
> <sup>3</sup>Department of Neurology, University of Pennsylvania, Philadelphia, PA 19104, USA

This repository contains analysis code and data related to research on characterizing individual differences in decision-making processes, specifically focusing on risk and benefit sensitivity in the context of cognitive enhancement attitudes.

## Related Manuscripts
- Fernandez, K., Erickson, B. A., Kable, J. W., Hamilton, R. H., & Medaglia, J. D. (2025). Evaluating Risk and Benefit Sensitivity for Cognitive Treatments. Journal of Cognitive Enhancement. https://doi.org/10.1007/s41465-025-00319-3
  
## Repository Contents
- Navigate to the `analysis` and `r_docs` directories to run the R scripts
- The `data` directory hosts raw and cleaned datasets. For details see [data dictionary](https://github.com/CogNeW/project_neuroethics_individual_difference/tree/master/data)
- The `images` directory contains figures 
- `output` houses intermediate data objects used in the analysis.

## Getting Started
To work with the datasets and analysis scripts:
1. Clone this repository:
   ```bash
   git clone https://github.com/CogNeW/project_neuroethics_individual_difference.git
   cd project_neuroethics_individual_difference
   ```

# Scripts of Interest (DDM Analysis)

## Time-Varying and standard DDM

The Time-Varying DDM was based off of [Chen et al.](https://www.nature.com/articles/s41598-022-15432-y).
The files are included under `analysis/timeDDM/Chen-2021/Scripts/Model Fitting`

1. `Neuroethics-stDDM.Rmd` - contains R code for runJags. Can run both the start-time DDM (stDDM / time-varying) and normal DDM on all participants as well as on individual participants.
2. `tssmHT_model_priors_M1.txt` - contains prior specifications for the stDDM. The specifications are hierarchical with both group-level and individual-level parameters being estimated.
3. `tssmHT_individual_model_priors_M1.txt` - contains non-hierarchical specifications for the stDDM. This file is used to computed parameters for individual participants. 
4. `tssmHT_modelWIENER_standardDDM.txt` - contains prior specificatons for standard 4-parameter DDM. The specifications are hierarchical with both group-level and individual-level parameters being estimated.
5. `tssmHT_modelWIENER_standardDDM_individual` - contains non-hierarchical specifications for the standard 4-parameter DDM. This file is used to computed parameters for individual participants. 

The other specification models are those included by Chen et al. in [their OSF](https://osf.io/k8g39/). These additional models were for additional assumptions and comparisons that were used in their paper. 

The results will be saved into the folder `analysis/timeDDM/Chen-2021/Scripts/Model Results`

## Regression Analysis

Our Regression analysis will be done in `analysis/DDM_STDDM_Regression.Rmd`. This analysis will read the stDDM results from `analysis/timeDDM/Chen-2021/Scripts/ModelResults/Neuroethics/20240313STDDMResults.csv` and perform regressions based on individual variables read from `data/Cleaned_Pilot/exclusions.2020.10.21.Rdata`.

# Description of File Directories and Raw Data

Included are:

r_docs: contains the functions, tidy scripts, RMarkdown, Stats analysis used throughout the project

analysis: RMarkdown files that constitute the final outputs used within the manuscript

data: contains the raw datasets prior to any filtering or measure score calculations. This will also have the cleaned versions of the data moving forward. 

images: contains the visualizations for subjects from a particular pilot cohort. To be used as reference for interpreting the findings

output: Intermediate data objects created in my analysis

ext: Any miscellaneous external files such as functions that take part in the analysis 

This raw data (raw.csv) contains the following columns:

- `user_id`	-- the unique subject i.d
- `rt`	-- recorded response time from each trial
- `key_press` -- the key that was pressed during the trial
- `trial_index` -- the total trial index. This includes instructional slides. 
- `time_elapsed`	 -- Total cumulative time elapsed within task
- `condition` --  the cognitive domain condition assigned to the subject
- `condition_per_trial`-- what cognitive domain condition the trial is. This should be the same for all trials within a given subject
- `total_trials` -- the total amount of trials a subject completed. This should be 75 for everyone. 
- `treatment_display`	-- which side the experimental treatment displayed on
- `risk` -- the number indicating how much risk is associated with the treatment for that trial
- `gain` -- the number indicating how likely you are to successfully treat the dysfunction using the treatment with that trial
- `no effect` -- the number indicating how likely it is that your treatment will cause no effect for treating the dysfunction
- `experimental_treatment_selected` -- the logical indicated if the subject selected the experimental treatment (TRUE) or no treatment (FASLE)
- `Progress` -- the percent completion of the survey battery
- `Duration..in.seconds.` -- the length of time taken to complete the survey in seconds
- `RecordedDate` -- the data the survey was completed
- `ResponseId` -- Recorded unique subject ID used to identify redundant subjects. That is, those who may have participated in the experiment twice
- `attention_check.X` -- the attention check items embedded. 1 = correct response, 0 = incorrect response
- `MFSQX` -- the MFSQX are the responses to the moral foundations survey (see preregistration). Each question is represented by the number following the identifier. e.g. MFSQ4 = question 4
- `BIG5_X` -- the BIG_X are the responses to the big 5 personality questionnaire (see preregistration). Each question is represented by the number following the identifier. e.g. BIG_4 = question 4
- `sex` -- biological sex: Male = 1, Female = 2, Other = 3.
- `age` -- age at time of survey
- `race_ethnicity` -- race or ethnicity of subject: Black = 1, White = 2, Hispanic or Latino = 3, Asian = 4, Native Hawaiian/Pacific Islander = 5, American Indian/Alaska Native = 6, Mixed = 7, Other = 8, Prefer not to say = 9.
- `education` -- highest level of education achieved: less than high school = 1, high school degree = 2,  bachelor = 3, masters = 4, Doctoral/Professional = 5
- `ReligionRaised` -- self-reported religion raised: Agnosticism, Athesism, Buddhism, Christianity, Hinduism, Islam, Judaism,Paganism, Other, and none.
- `ReligionNow` -- self-reported current religion: Agnosticism, Athesism, Buddhism, Christianity, Hinduism, Islam, Judaism,Paganism, Other, and none.
- `Religiosity` -- recoded likert scale of how religious the participant reported ranging from Extremely Religious to Not Religious
- `Spirituality` -- recoded likert scale of how spiritual the participant reported ranging from Extremely Religious to Not Religious
- `AttendChurch` -- How often the participant reported attending church: Never, Once a year or less, A few times a year, A few times a month, Once a week, or More than once a week. 
- `ReligiousActivity` -- How often the participant reported engaging in religious activities: Never, Once a year or less, A few times a year, A few times a month, Once a week, or More than once a week.
- `meditation` -- How often the participant reported meditating: Never, Once a year or less, A few times a year, A few times a month, Once a week, or More than once a week.
- `RelationshipStatus` -- relationship status of participant: Single, In a Relationship, Married, Divorced, or Widowed.
- `Income` -- Reported income in US dollars rounded to the nearest dollar.
- `Major` -- if applicable, reported college major.
- `YearOfEducation` -- number of total years of education.
- `EmploymentStatus` -- if applicable, reported current occupation.
- `Political_Party` -- reported political party affiliation: Democrat, Republican, Independent, or Other. 
- `HumanMind` -- If you had to choose, do you think that human minds are more like computers or souls?
- `Optimize` -- Y/N: Have you ever used a device to optimize cognition? 
- `republican` -- the degree to which republican describe your political identity from weak to strong on a 7 point likert scale. 
- `democrat` -- the degree to which democrat describe your political identity from weak to strong on a 7 point likert scale.
- `independent` -- the degree to which independent describe your political identity from weak to strong on a 7 point likert scale.
- `conservative` -- the degree to which conservative describe your political identity from weak to strong on a 7 point likert scale.
- `liberal` -- the degree to which liberal describe your political identity from weak to strong on a 7 point likert scale.
- `POWER` -- the importance of the power as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `SELF.DIRECTION` -- the importance of the self direction as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `ACHIEVEMENT` -- the importance of the achievement as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `HEDONISM` -- the importance of the hedonism as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `STIMULATION` -- the importance of the stimulation as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `UNIVERSALISM` -- the importance of the universalism as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `BENEVOLENCE` -- the importance of the benevolence as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `TRADITION` -- the importance of the tradition as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles,  1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `CONFORMITY` -- the importance of the conformity as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles, 1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `SECURITY` -- the importance of the security as a life-guiding principle using an 8-point scale in which 0 indicates that the value is opposed to your principles, 1 indicates that the value is not important for you, 4 indicates that the values is important, and 8 indicates that the value is of supreme importance for you. 
- `MC_X` -- the MC_X are the responses to the moral commitments task (see preregistration). Each question is represented by the ethical area under consideration. e.g. MC_Lying = moral commitments question concerning what the participant thought was the best representation of what is right or wrong to do.
- `FSX` -- Interpersonal reactivity index (IRI) Fantasy subscale. A 5-point Likert scale ranging from “Does not describe me well” to “Describes me very well” (for more information see preregistration).
- `ECX` -- Interpersonal reactivity index (IRI) Empathic Concern subscale. A 5-point Likert scale ranging from “Does not describe me well” to “Describes me very well” (for more information see preregistration).
- `PTX` -- Interpersonal reactivity index (IRI) Perspective Taking subscale. A 5-point Likert scale ranging from “Does not describe me well” to “Describes me very well”(for more information see preregistration).
- `PDX` -- Interpersonal reactivity index (IRI) Personal Distress subscale. A 5-point Likert scale ranging from “Does not describe me well” to “Describes me very well”(for more information see preregistration).
- `SNSX` --  Subjective Numeracy Scale (SNS) 8-item scale question (for more information see preregistration).
- `BerlinX` -- Berlin Numeracy task (for more information see preregistration).
- `GRiPS_X` -- self-report measure of general risk propensity. 8 items using a 5-point Likert scale 1 = strongly disagree to 5 = strongly agree. (for more information see preregistration).
- `Q1802_X` -- self reported data quality check questions. scaled 1 = Disagree Strongly to 7 = Agree Strongly (for more information see materials). 


The cleaned data frame (clean_data_long.csv) is a long format data frame containing the following columns:

- `user_id` --  Unique subject identifier
- `risk` -- the number indicating how much risk is associated with the treatment for that trial
- `gain` -- the number indicating how likely you are to successfully treat the dysfunction using the treatment with that trial
- `no effect` -- the number indicating how likely it is that your treatment will cause no effect for treating the dysfunction

   
## Contact
For any inquiries, please contact the corresponding author, John D. Medaglia, at [jdm383@drexel.edu](mailto:jdm383@drexel.edu).
