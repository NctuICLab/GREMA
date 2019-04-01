# GREMA: Modelling of emulated gene regulatory networks with confidence levels based on evolutionary intelligence to cope with the underdetermined problem
This work proposes an evolutionary modelling algorithm (EMA) that is based on evolutionary intelligence, including both crowd wisdom and an evolutionary strategy, to cope with the underdetermined problem. EMA uses an intelligent genetic algorithm to solve the large-scale parameter optimisation problem. 
An EMA-based method, GREMA, infers a novel type of gene regulatory network with confidence levels for every inferred regulation. The higher the confidence level, the more accurate the inferred regulation. GREMA gradually determines the regulations of an eGRN with confidence levels in descending order using either an S-system or a Hill function-based ordinary differential equation model. 

Setup and Data Format
============================
There are 2 input tsv files: time-series data and domain knowledge data.
The format of time-series data has 2 parts:
 - repeat_number= (ex:repeat_number=3)
 - time-series profile (rep1 <tab> gene_name <tab> gene_expression)
  
