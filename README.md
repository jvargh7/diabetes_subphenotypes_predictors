# subphenotypes_predictors
Predictors of subphenotypes in cohort studies and EHRs


 - Study Objective: 
 To study the associations of key pathophysiological biomarkers with membership in data-driven type 2 diabetes subtypes, 
 compared to individuals who did not develop type 2 diabetes.
 
 
 - Study population:
 Newly diagnosed T2D + No T2D participants from four US cohorts (CARDIA, DPP/OS, JHS, MESA)
 
 
 - Repo structure
 
 i. `preprocessing`:
 
 `dsppre01a-h` Keep new T2D cases + No T2D in each cohort, recode variables.
 `dsppre01` Harmonized all cohorts, identify new T2D cases, calculate HOMA2 indices.
 
 
 ii. `hyperc3`:
 
 Multiple imputation (by cohort) code run from Emory HPC cluster.
 
 
 iii. `analysis`:
 
 `dspan01` Restrict the analytic sample to those with first and last visits with available HbA1c and BMI. Each selected participant is required to have at least two visits.
 `dspan02` Explore the trajectories of four biomarkers in the 15y preceding T2D diagnosis (Figure 2).
 `dspan03` TDCM survival analysis to study the associations between biomarkers and subtypes (Figure 1).
 `dspan04` Descriptive analysis (Table 1).
 `dspan05` Cox PH, SDH, multinomial models using baseline data (sensitivity analysis).
 `dspan06` Evaluate Cox PH model discrimination and calibration (sensitivity, specificity, F1 score, C index; calibration slope) using baseline data.
 `dspan07` Evaluate Cox PH model calibration (calibration slope) for each imputed baseline datasets.
 
 iv. `sensitivity analysis`:
 
 `dspse01a-b` Apply IPW to address potential selection bias from missing HbA1c and BMI at baseline or follow-up.
 `dspse02a-b` Impute datasets separately by cohort, pool results.
 `dspse03`    Remove DPP intervention arm.
 `dspse04`    Add HDL and TG.
 `dspse05`    Add other biomarkers to the trajectory plot.
 

 
 