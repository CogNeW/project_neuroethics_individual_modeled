# Code for 'The Dynamics of Neuroethical Decisions Depend on Psychological Traits'

> **Authors:**  
> Brian Kim<sup>1</sup>, Kiante Fernandez<sup>1,2</sup>, John D. Medaglia<sup>1,3</sup>  
> <sup>1</sup>Department of Applied Cognitive and Brain Sciences, Drexel University, Philadelphia, PA, 19104, USA  
> <sup>2</sup>Department of Psychology, University of California, Los Angeles, CA 90095, USA  
> <sup>3</sup>Department of Neurology, University of Pennsylvania, Philadelphia, PA 19104, USA

This repository contains analysis code and data related to research on characterizing individual differences in decision-making processes, specifically focusing on risk and benefit sensitivity in the context of cognitive enhancement attitudes.

## Terminology

The manuscript reports these constructs as **harm** and **cure** probabilities. The
code and stored data predate that change and still use the earlier *risk* / *benefit*
(or *gain*) vocabulary internally. The mapping is one-to-one:

| Manuscript term | Code / data term | Variable |
| --- | --- | --- |
| Harm probability, harm sensitivity | Risk probability, risk sensitivity | `b_risk`, `risk` |
| Cure probability, cure sensitivity | Benefit / gain probability, gain sensitivity | `b_gain`, `gain` |
| Cure first / harm first (latency) | Benefits first / risks first | `time` |

Only the strings that appear in rendered figures were updated to the new terms.
Variable names, column names, and `.RData` object names were deliberately left
unchanged so that saved model fits and joins continue to work. Note that
`Risk Taking (GRiPS)` is the name of a published scale and is not part of this
mapping.

## Related Manuscripts
- Fernandez, K., Erickson, B. A., Kable, J. W., Hamilton, R. H., & Medaglia, J. D. (2025). Evaluating Risk and Benefit Sensitivity for Cognitive Treatments. Journal of Cognitive Enhancement. https://doi.org/10.1007/s41465-025-00319-3
  
## Repository Contents
- The `analysis` folder hosts our main analysis and figure generation scripts
   - `finalDDM/` contains the finalized versions used for the paper
   - `archive/` contains intermediate draft versions used for development
- The `data` directory hosts raw and cleaned datasets.

# Scripts of Interest (DDM Analysis)

## Time-Varying and standard DDM

The Time-Varying DDM was based off of [Chen et al.](https://www.nature.com/articles/s41598-022-15432-y).

The main analysis script is included under `analysis/finalDDM/kim_analysis_stddm`. The file references two jags model configurations, specifically `analysis/finalDDM/stddm_model_M0.txt` and `analysis/finalDDM/stddm_model_M1.txt`, which specify the simple and time-varying (start-time) DDM respectively.




## Regression Analysis

Our Regression analysis will be done in `analysis/DDM_STDDM_Regression.Rmd`. This analysis will read the stDDM results from `analysis/timeDDM/Chen-2021/Scripts/ModelResults/Neuroethics/20240313STDDMResults.csv` and perform regressions based on individual variables read from `data/Cleaned_Pilot/exclusions.2020.10.21.Rdata`.

## Participant Exclusions and Robustness (reviewer response)

`analysis/finalDDM/exclusion_robustness.R` addresses the reviewer question about
the 240 participants excluded from the recruited sample of 616. It is standalone:
it reads `data/raw_mturk_data/raw.csv`,
`data/Cleaned_Pilot/exclusions.2020.10.21.RData`,
`data/FINAL_RT_EXCLUDED_DATA_20241014.RData`,
`analysis/finalDDM/FINAL_PARTICIPANT_LIST.RData` and
`analysis/finalDDM/final_ddm_ideological_combined.RData`, and writes every table,
figure and a full text transcript to
`analysis/finalDDM/Figures/ExclusionRobustness/`. Run it with
`Rscript analysis/finalDDM/exclusion_robustness.R` (needs `dplyr`, `tidyr`,
`ggplot2`, `psych`, `arm`); it takes a couple of minutes.

What it produces:

1. **Exclusion cascade** (616 -> 610 -> 425 -> 382 -> 376), reconciling exactly to
   the n = 240 reported in the manuscript. The 240 is a total across several
   criteria: 6 with no usable record, 185 removed by the prior-study criteria
   applied upstream (Fernandez et al., 2025), 43 for choice profile, and 6
   dropped between the DDM fit and the analysis sample.
2. **The response criterion as implemented.** The Methods text says participants
   without both response types were excluded. The rule actually applied is a
   10%/90% band on the Recommend rate *computed after RT trimming*: keep
   participants with `0.10 < p < 0.90`. The script verifies this reproduces the
   382/43 split with zero misclassifications. Two consequences are reported and
   both should go in the Methods:
   - it is a proportion rule, not a count rule, so minority-response *counts*
     overlap between kept and dropped (4 vs 7) even though proportions do not;
   - because the proportion is computed post-trimming, inclusion is jointly
     determined by the response profile and by which trials the RT filter
     removed. Re-applying the same band to untrimmed trials flips 10 of 425
     participants (2.4%). That interaction is worth disclosing.

   It is *not* the case that trimming stripped minority responses until people
   became degenerate: trimming removes a median of one minority response among
   the dropped, none were strictly one-sided beforehand, and only 2 of 43 are
   afterwards. The cut is a deliberate near-uniformity threshold.
3. **Choice-profile breakdown** of every excluded participant, classified on the
   same post-trimming rate the criterion used, so the "strong preferences"
   hypothesis can be checked directly. Within the choice-profile tranche the
   split is 21 near-always-reject vs 22 near-always-recommend.
4. **Trait comparison**, analyzed vs excluded, on the measures whose scoring keys
   reproduce the published composites exactly (GRiPS, Berlin numeracy, subjective
   numeracy, BFI-2 domains, religiosity, spirituality, age): Welch tests, Hedges'
   g with CIs, BH-FDR, Kolmogorov-Smirnov, and TOST equivalence tests against
   d = 0.35. Also run separately for the choice-profile tranche and for
   reject-side vs recommend-side within the excluded.
5. **Response-style diagnostics** (acquiescence, within-person response
   variability, Berlin accuracy, survey duration) that distinguish a substantive
   trait difference from low-effort questionnaire responding.
6. **Model-free robustness check.** Per-participant penalized logistic
   regressions of choice on harm and cure probability, which stay estimable for
   near-degenerate responders where the stDDM does not, letting the key trait
   associations be re-run with the excluded participants added back in.
7. **Threshold sensitivity.** The same regression bank re-fit at bands from
   0.00/1.00 (keep everyone) through 0.20/0.80 (twice as strict as published),
   showing whether the trait associations depend on where the cut was drawn.

Note that section 4b of the script flags a scoring issue that is independent of
the exclusion question: the published BFI-2 domain scores are plain item means
taken without reverse-scoring, which leaves them dominated by acquiescence (mean
inter-domain r = .80). The exclusion analysis reports results under both
scorings, but the Big Five results in the manuscript rest on the unreversed
scores and should be re-derived against the upstream cleaning pipeline.

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

   
## BFI-2 Scoring Correction

The five BFI-2 domain columns (`Extraversion`, `Agreeableness`, `Conscientiousness`,
`Neuroticism`, `Openness`) were originally stored as plain means over all twelve items
of each domain, **with no reverse-scoring**. Each BFI-2 domain is balanced across
forward- and reverse-worded items, so the reverse-worded items cancelled the
forward-worded ones and what remained was largely acquiescent responding.

Evidence, all reproducible from files in this repository:

- Items in `data/raw_mturk_data/raw.csv` are raw 1-5 responses. Opposite-worded items
  correlate negatively, e.g. r(`BIG5_1` "Is outgoing, sociable", `BIG5_16` "Tends to be
  quiet") = -.42, while r(`BIG5_1`, `BIG5_46` "Is talkative") = +.67.
- The stored columns reproduced `rowMeans()` over the unreversed items to machine
  precision for all 382 analyzed participants.

| | as stored | corrected |
| --- | --- | --- |
| Mean inter-domain r | +.68 | +.02 |
| Cronbach's alpha | .62 - .70 | .79 - .88 |
| Mean domain SD | 0.45 | 0.80 |

`analysis/finalDDM/fix_bfi2_composites.R` recomputes the five columns in
`data/FINAL_RT_EXCLUDED_DATA_20241014.RData` using the published key (Soto & John,
2017; six of the twelve items reverse-worded in every domain). It asserts the
pre-fix state before writing and copies the original to
`data/FINAL_RT_EXCLUDED_DATA_20241014.PRE_BFI2_FIX.RData`. Note that `.RData` files
are excluded by `.gitignore`, so this correction must be re-run against any local or
archived copy of that file.

### This does not change any manuscript result

The composites are carried into `final_table` but never enter a reported analysis;
they appear in neither `independent_vars` set in
`analysis/finalDDM/DDM_STDDM_Regression.Rmd`. Every published Big Five result is
item-level ridge importance -- `abs(coef(ridge_model)[-1])` with `standardize = TRUE`
(line 536) -- and reverse-scoring an item flips the sign of its coefficient while
leaving the magnitude unchanged. Verified directly on these items: maximum difference
in importance between unreversed and correctly keyed inputs is 4e-16, with all 30
reverse-keyed items flipping sign. The Agreeableness / `time` permutation result
(FDR p = .039) is therefore unaffected.

### Caveat for any signed or directional reuse

Importance is sign-free, but **domain-averaged signed coefficients are not stable
under this correction**. `melted_coefs` in `RidgeResults.RData` stores signed ridge
beta weights; averaging them within a Big Five domain averages six correctly signed
items against six inverted ones. Across the five domains and the `time` / `b_gain` /
`b_risk` parameters, 7 of 15 such domain means change sign once the keying is
corrected, and those that do not change sign can change substantially in magnitude
(Agreeableness / `time`: -7.1e-5 as stored vs -9.3e-4 corrected). `melted_coefs` is
written but never read by any script in this repository. Any directional claim about
a Big Five domain should be re-derived from correctly keyed items rather than taken
from a signed average over the stored ones. Objectively scored measures (Berlin
Numeracy, scored correct/incorrect) are not affected.

### Composite columns in the pilot CSVs

The same unreversed composites also appear in five tracked pilot files:
`data/Cleaned_Pilot/clean_data_long.csv`, `clean_data_longV2.csv`, `rawbeta.csv`,
`withoutselfreportExlcusion.csv`, and `meancentered.csv` (mean inter-domain r = +.42
to +.72). These are pilot data and are not corrected by the script above.

## Contact
For any inquiries, please contact the corresponding author, John D. Medaglia, at [jdm383@drexel.edu](mailto:jdm383@drexel.edu).
