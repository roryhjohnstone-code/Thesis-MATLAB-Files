# Fibre Bragg Grating (FBG) Based Structural Health Monitoring for Wind Turbine Blades

This repository contains the MATLAB simulation framework developed for the MSc thesis  
**"Fibre Bragg Grating (FBG) Based Structural Health Monitoring for Wind Turbine Blades"**  
(University of Aberdeen, 2025).

The system models Fibre Bragg Grating (FBG) optical sensors embedded along a turbine blade to detect structural faults such as delamination, cracks, impact, and sensor drift.  
It combines adaptive thresholding, Kalman filtering, and machine-learning classification inside an interactive MATLAB GUI, with an additional benchmarking mode for automated testing.

---

Quick Start

Open MATLAB and run one of the following scripts from the **`demos`** folder:

1. `run_demo` — Launch the main GUI
```matlab
cd demos
run_demo
```
This opens the interactive FBG SHM GUI (fbg_shm_gui.m), which provides:
Real-time plots of strain and fault detection for all sensors
Controls to inject faults (Drift / Crack / Delamination / Impact)
Kalman filtering and ML classifier toggles
Adaptive thresholding with hysteresis
3-D animated blade visualisation
Real-time diagnostic table and export functions

 \run_benchmark_suite - executes the SHM system automatically, without the GUI.
```matlab
cd demos
run_benchmark_suite
```
It performs a full benchmarking sweep by:
Testing four fault types: Drift, Delamination, Crack, Impact
Varying signal-to-noise ratio (SNR) from 50 dB to –5 dB
Enabling/disabling Machine-Learning and Kalman Filtering options
Repeating each configuration several times for statistical reliability
For each case, it logs:
Detection accuracy per sensor and overall
Precision, recall, and F1-score
Median detection delay (seconds) 
All results are automatically written to:
```matlab
demos/exports/benchmark_YYYYMMDD_HHMMSS.csv
```

Machine Learning The machine learning subsystem enables automated fault classification, It provides a complete pipeline for generating training data, extracting features, and applying supervised classification
**`train_and_cache_model.m`** 
  Trains a supervised classifier — a *bagged ensemble of decision trees* — using statistical features extracted from simulated strain data. 
  The trained model is cached locally (`bestClassifier.mat`) for reuse, ensuring reproducibility and fast execution across simulation runs.

- **`applyMLPerSensor.m`** 
  Applies the trained model to each sensor's data stream, returning per-sensor fault predictions (`Normal`, `Drift`, `Delamination`, `Crack`, `Impact`) and confidence scores. 
  This allows multi-sensor fusion of classification results and enables visualization of fault propagation across the turbine blade.

- **`generateTrainingData.m`** 
  Generates synthetic, labelled datasets representing both *healthy* and *faulty* structural conditions. 
  Each class includes multiple noise levels and fault amplitudes to improve model robustness and generalisation.

- **`extractFBGFeatures.m`** 
  Extracts key time domain and statistical features — such as mean, variance, RMS, kurtosis, skewness, and peak-to-peak values from each sensor's strain signal. 
  Ensures consistent feature alignment between training and runtime classification.
Together, these scripts form a modular pipeline for **data generation → feature extraction → model training → real-time classification**, enabling the MATLAB GUI to perform automated fault detection and diagnosis using a reproducible bagged-tree ensemble model.

Repository Strucutre
src/       → MATLAB source code (GUI, simulation, ML, Kalman, visualisation)
demos/     → Runnable scripts (run_demo.m, run_benchmark_suite.m)
MATLAB_TOOLBOXES.md → Required MATLAB toolboxes

Requirements
MATLAB R2024b or newer
Signal Processing Toolbox
Statistics and Machine Learning Toolbox
Communications Toolbox

Licensce
MIT License © Rory Johnstone 2025