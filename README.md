Datasets and R codes associated with the manuscript titled "Predator learning can resolve the paradox of local warning signal diversity"
Below is a list of files:

1. attack_prob_curated.csv: dataset used for the meta-analysis on the variation in final attack probability  
**author:** first author of the publication  
**journal:** the journal in which the study was published  
**volume:** the volume of the journal  
**species:** study species  
**prey_type:** type of prey (almond, mealworm, chicken feed, etc)    
**prey_presentation:** the manner in which prey were presented to the predators (sequential or simultaneous)  
**cs:** type of conditioned stimulus, or the ways in which profitable and control prey consistently differed  
**cs2:** same as cs, coded as either color, pattern, color+pattern, or color+other (size or odor)  
**cs_property:** additional comments or notes for cs  
**us:** cause of unprofitability for artificial prey, recorded as either the type of bitter chemicals (when applicable) or prey defense
**us2:** same as us, but coded as either bitterness (when prey unprofitability was from the use of bitter chemicals) or toxin (when prey extracts were used)  
**us_concentration:** the original concentration of bitter chemicals used (expressed as weight percentage)  
**us_concentration_cal:** calibrated bitterness, accounting for the fact that denatonium benzoate was 20 times more bitter than quinine derivatives  
**attack_prob:** final attack probabilities of predators against unprofitable prey  
**way_of_learning:** how predators learn, recorded as either direct or indirect (i.e., through social learning or generalization)  
**prey_presentation2:** same information as prey_presentation, but coded as either 1 or 2 for statistical purposes  
**species2:** same as species, but coded as 1-4 for statistical purposes  
**cs3:** same as cs2, but coded as 1-4 for statistical purposes    
**us3:** same as us2, but coded as either 1 or 2 for statistical purposes  
  
2. P&F.csv: Final attack probabilities and forgetting rates from selected predators  
**species:** name of the predator species  
**taxa:** taxomomic group of the predator species  
**P:** final attack probabilities  
**F:** forgetting rates  
**source:** the study from which I obtained the data  
**note:** miscellaneous notes  
  
3. learning_IMB.R: R code for performing individual-based simulations of predator-prey interactions

4. learning_ODE.R: R code for performing population-level simulations of predator-prey interactions  
  
5. learning_functions.R: containing custom functions necessary for performing simulations  
  
6. meta-analysis.R: R code for statistical meta-analysis  
  
7. HighstatLibV10.R: R code from Zuur et al. (2009) containing custom functions for GL(M)M model validation

8. Supplementary Materials_ver2.docx: Supplementary tables, figures,and references accompanying the main paper.
