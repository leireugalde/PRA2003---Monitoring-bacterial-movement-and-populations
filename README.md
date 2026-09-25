# PRA2003---Monitoring-bacterial-movement-and-populations
Leire Ugalde Salvá - i6384470

## Overview ##
This project analyzes bacterial data tracking to study the movement and proliferation of different bacteria under specific nutrient or stress conditions. It analyzes changes in behavior in bacterial strains and their genetic variants. For each event, the data is read under and used to answer the three questions presented below. 

## Answer the following questions with the given dataset on bacterial movement and populations ##
1. What are the average counts of each bacterial stain and their statistical uncertainties?
2. Is there any asymmetry between the normal and the mutant strain?
3. Is there any asymmetry as a function of their momentum?

## Files ##
• (PRA2003Week2.R) reads one event from the data file (output-Set0.txt) and calculates the momentum magnitude of each bacteria


• (PRA2003Week3.R) reads the average bacteria per event together with the statistical uncertainty (default output-Set1.txt, but allows user input) ---- This was the first attempt at the code, it had errors so it is not used for following weeks

• (PRA2003Week3Corrected.R) correctly reads average bacteria per event, together with statstical uncertainty (default output-Set1.txt). This code is used in the following weeks

• Following files will be defined in upcoming weeks

## Data format ##
• Data file (output-Set0.txt) is structured as 
    -Header line: event_id  n_particles
    -One line per particle: px py pz code 
        --> Where px/py/pz are momentum components 
        --> Where code identifies the strain (normal or mutant)

• The rest of data is formated in the same manner (output-Set1.txt, output-Set2.txt...)

## Results ##
To be filled in once complete analysis is completed 


