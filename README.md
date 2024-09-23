# TD-MRC
Bivariate **T**emporal **D**ependence via **M**ixtures of **R**otated **C**opulas

This repository contains various models and results related to mixture models and Clayton copulas. Below is a summary of the contents of each folder and its purpose.

## Data
  - This folder contains the real data. Full data and split data for cross-validation.

## Implementations

### 1. **MixClayton**
   - Includes files implementing the mixture Clayton model.

### 9. **SingleClayton**
   - Includes implementation and results for the single Clayton model.
     
### 4. **Gaussian**
   - This folder holds files and scripts related to the Gaussian models.

## Results

### 1. **Simulation_Result**
   - This folder contains the results from the three different models. Simulated data is from a dynamic mixture of Clayton.

### 2. **Fitting_Result**
   - Contains the results from fitting models to real data. The files provide full insights into how the models fit under different dependence structures and strengths.
     
### 3. **CV_Result**
   - Contains cross-validation results from the applied models. These results are used to evaluate the performance of different models based on their predictive accuracy.

### 4. **Prediction_Result**
   - Contains prediction results from the applied models.  Mainly using LPS (log predictive score) to evaluate Their performance were evaluated mainly by LPS (log predictive score).

### 5. **Result**
   - Contains some old results which can be ignored.
