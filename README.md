# Environmental Risk Score (ERS)
This github repository consists of R code for construction of the "Environmental Risk Score (ERS)"-- a tool designed for examining the risk of exposure to multi-pollutants in epidemiologic research, in the following papers:

## Environmental Risk Score as a New Tool to Examine Multi-Pollutants in Epidemiologic Research: An Example from the NHANES Study Using Serum Lipid Levels
**ABSTRACT**

OBJECTIVE: A growing body of evidence suggests that environmental pollutants, such as heavy metals, persistent organic pollutants and plasticizers play an important role in the development of chronic diseases. Most epidemiologic studies have examined environmental pollutants individually, but in real life, we are exposed to multi-pollutants and pollution mixtures, not single pollutants. Although multi-pollutant approaches have been recognized recently, challenges exist such as how to estimate the risk of adverse health responses from multi-pollutants. We propose an “Environmental Risk Score (ERS)” as a new simple tool to examine the risk of exposure to multi-pollutants in epidemiologic research.

METHODS and RESULTS: We examined 134 environmental pollutants in relation to serum lipids (total cholesterol, high-density lipoprotein cholesterol (HDL), low-density lipoprotein cholesterol (LDL) and triglycerides) using data from the National Health and Nutrition Examination Survey between 1999 and 2006. Using a two-stage approach, stage-1 for discovery (n = 10818) and stage-2 for validation (n = 4615), we identified 13 associated pollutants for total cholesterol, 9 for HDL, 5 for LDL and 27 for triglycerides with adjustment for sociodemographic factors, body mass index and serum nutrient levels. Using the regression coefficients (weights) from joint analyses of the combined data and exposure concentrations, ERS were computed as a weighted sum of the pollutant levels. We computed ERS for multiple lipid outcomes examined individually (single-phenotype approach) or together (multi-phenotype approach). Although the contributions of ERS to overall risk predictions for lipid outcomes were modest, we found relatively stronger associations between ERS and lipid outcomes than with individual pollutants. The magnitudes of the observed associations for ERS were comparable to or stronger than those for socio-demographic factors or BMI.

CONCLUSIONS: This study suggests ERS is a promising tool for characterizing disease risk from multi-pollutant exposures. This new approach supports the need for moving from a single-pollutant to a multi-pollutant framework.

##### **CITATION**:
Park, S.K., Tao, Y., Meeker, J.D., Harlow, S.D. and Mukherjee, B., 2014. Environmental risk score as a new tool to examine multi-pollutants in epidemiologic research: an example from the NHANES study using serum lipid levels. PloS one, 9(6), p.e98632.
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098632

## Construction of environmental risk score beyond standard linear models using machine learning methods: application to metal mixtures, oxidative stress and cardiovascular disease in NHANES
**ABSTRACT**

BAKCGROUND: There is growing concern of health effects of exposure to pollutant mixtures. We initially proposed an Environmental Risk Score (ERS) as a summary measure to examine the risk of exposure to multi-pollutants in epidemiologic research considering only pollutant main effects. We expand the ERS by consideration of pollutant-pollutant interactions using modern machine learning methods. We illustrate the multi-pollutant approaches to predicting a marker of oxidative stress (gamma-glutamyl transferase (GGT)), a common disease pathway linking environmental exposure and numerous health endpoints.

METHODS: We examined 20 metal biomarkers measured in urine or whole blood from 6 cycles of the National Health and Nutrition Examination Survey (NHANES 2003–2004 to 2013–2014, n = 9664). We randomly split the data evenly into training and testing sets and constructed ERS’s of metal mixtures for GGT using adaptive elastic-net with main effects and pairwise interactions (AENET-I), Bayesian additive regression tree (BART), Bayesian kernel machine regression (BKMR), and Super Learner in the training set and evaluated their performances in the testing set. We also evaluated the associations between GGT-ERS and cardiovascular endpoints.

RESULTS: ERS based on AENET-I performed better than other approaches in terms of prediction errors in the testing set. Important metals identified in relation to GGT include cadmium (urine), dimethylarsonic acid, monomethylarsonic acid, cobalt, and barium. All ERS’s showed significant associations with systolic and diastolic blood pressure and hypertension. For hypertension, one SD increase in each ERS from AENET-I, BART and SuperLearner were associated with odds ratios of 1.26 (95% CI, 1.15, 1.38), 1.17 (1.09, 1.25), and 1.30 (1.20, 1.40), respectively. ERS’s showed non-significant positive associations with mortality outcomes.

CONCLUSIONS: ERS is a useful tool for characterizing cumulative risk from pollutant mixtures, with accounting for statistical challenges such as high degrees of correlations and pollutant-pollutant interactions. ERS constructed for an intermediate marker like GGT is predictive of related disease endpoints.

##### **CITATION**:
Park, S.K., Zhao, Z. and Mukherjee, B., 2017. Construction of environmental risk score beyond standard linear models using machine learning methods: application to metal mixtures, oxidative stress and cardiovascular disease in NHANES. Environmental Health, 16(1), p.102.
https://ehjournal.biomedcentral.com/articles/10.1186/s12940-017-0310-9

##### **SAMPLE R CODE**: ERS_Environmental Health_2017_.R


## Associations of cumulative exposure to heavy metal mixtures with obesity and its comorbidities among U.S. adults in NHANES 2003–2014
**ABSTRACT**

BACKGROUND: Some heavy metals (e.g., arsenic, cadmium, lead, mercury) have been associated with obesity and obesity comorbidities. The analytical approach for those associations has typically focused on individual metals. There is a growing interest in evaluating the health effects of cumulative exposure to metal mixtures.

OBJECTIVES: We utilized our Environmental Risk Score (ERS), a summary measure to examine the risk of exposure to multi-pollutants in epidemiologic research, to evaluate the associations of cumulative exposure to a mixture of correlated heavy metals with obesity and its comorbidities including hypertension, and type-2 diabetes mellitus (T2DM) while accounting for high degree correlations and interactions among metal mixtures components.

METHODS: We examined blood and urinary markers of 18 heavy metals among 9537 adults in NHANES 2003–2014. We randomly split data into a training set for the construction of ERS (n = 6675) and a testing set for the evaluation of its statistical performance (n = 2862). ERS of heavy metal mixtures was computed for waist circumference using adaptive elastic-net (AENET) with 189 predictors including 18 main effects, 18 squared terms, and 153 pairwise interactions of heavy metals. Regression analyses with complex survey designs were performed to assess the associations of ERS with other obesity measures, hypertension and T2DM.

RESULTS: 7 main effects (blood lead, blood cadmium, blood mercury, and urinary markers of monomethylarsonic acid (MMA), barium, mercury and thallium), 4 squared terms (blood cadmium, urinary cadmium, urinary antimony and urinary tungsten), and 7 pairwise interactions (blood lead & urinary cadmium, blood lead & urinary MMA, blood lead & urinary uranium, urinary cadmium & urinary MMA, urinary dimethylarsinic acid (DMA) & urinary tungsten, urinary MMA & urinary cobalt, and urinary lead & urinary antimony) were selected by AENET for construction of ERS of waist circumference-related metal mixtures. An increase in ERS from 10th percentile to 90th percentile in the overall study population was significantly associated with 4.50 kg/m2 (95% CI: 4.06, 4.94) higher BMI, 4.16 mm (95% CI: 3.56, 4.76) higher skinfold thickness, and 4.11 kg (95% CI: 0.83, 7.40) higher total body fat, independent of age, sex, race/ethnicity, education, smoking status, physical activity and NHANES cycle (Ps < 0.05). Significant associations of ERS with both hypertension and T2DM were also observed (Ps < 0.05).

CONCLUSIONS: Our study suggests that cumulative exposure to heavy metals as mixtures is associated with obesity and its related chronic conditions such as hypertension and T2DM. Additional research is needed to confirm these findings in longitudinal settings.

##### **CITATION**:
Wang, X., Mukherjee, B. and Park, S.K., 2018. Associations of cumulative exposure to heavy metal mixtures with obesity and its comorbidities among US adults in NHANES 2003–2014. Environment international, 121, pp.683-694.
https://www.sciencedirect.com/science/article/pii/S0160412018312650

##### **SAMPLE R CODE**: ERS_Obesity_NHANES_Environment_International_2018.R 
