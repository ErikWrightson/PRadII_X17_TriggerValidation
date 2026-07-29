# Trigger Validation Analysis Script Framework

## Overview

This repository contains analysis tools for processing, validating, and evaluating trigger performance using reconstructed detector data. The workflow is designed to process large numbers of ROOT files, combine trigger validation outputs, and generate diagnostic plots that characterize trigger efficiency, timing behavior, and spatial response across HyCal.

The analysis chain consists of:

1. Running the trigger validation analysis over reconstructed data files.
2. Combining individual ROOT outputs into a single merged file.
3. Producing trigger performance summary plots.
4. Evaluating trigger efficiency as a function of position across the HyCal calorimeter.

---

# Workflow Overview

```text
Raw/Reconstructed Data
          |
          v
+-----------------------+
| Run_triggerStudy1.py  |
| Batch Trigger Analysis|
+-----------------------+
          |
          v
TrigVal_<RUN>_All.root
          |
          +-----------------------------+
          |                             |
          v                             v
+-------------------+        +----------------------+
| Trigger Summary   |        | XY Trigger Efficiency|
| Plotting Script   |        | Analysis Script      |
+-------------------+        +----------------------+
          |
          v
 Trigger Performance PDFs
```

---

# General Notes

* The plotting scripts operate on ROOT files produced by the `TrigVal` trigger validation analysis.
* Histogram naming conventions must remain consistent for automatic histogram discovery.
* Trigger efficiency calculations use binomial statistical uncertainties.
* The XY efficiency analysis evaluates detector-position-dependent trigger response across HyCal.
* The full workflow provides both global trigger validation and spatial performance characterization.

This README serves as the central reference for the trigger validation analysis framework.

Below is a list of a separate README section for each script that can be useed in the analysis workflow of the X17 Clustering trigger analysis versus for the random trigger runs.

1. [Batch Runner Script](#1-batch-trigger-validation-runner)
2. [Summary Plotting for Merged Validation Histograms](#2-trigger-validation-summary-plotter)
3. [HyCal XY Trigger Efficiency Analysis and Plotting](#3-hyCal-xY-trigger-efficiency-analysis)
---
# 1. Batch Trigger Validation Runner

## Overview

`Run_triggerStudy1.py` is a Python utility that automates the execution of the `TrigVal` trigger validation analysis over all reconstructed data files for a given run. The script processes files in parallel, merges the resulting ROOT files into a single output, and removes the intermediate per-file ROOT files upon successful completion.

This utility is intended to simplify large-scale trigger validation analyses while reducing total processing time through parallel execution.

## Workflow

Given a run number, the script performs the following steps:

1. Searches the input data directory for all reconstructed ROOT files corresponding to the specified run.
2. Creates the output directory if it does not already exist.
3. Launches one `TrigVal` analysis job for each input file.
4. Executes multiple jobs simultaneously (up to a configurable maximum number of parallel processes).
5. Waits for all jobs to finish.
6. Merges the individual output ROOT files into a single ROOT file using `hadd`.
7. Deletes the intermediate per-file ROOT files after a successful merge.
8. Prints a summary of successful and failed jobs.

## Input

The script expects reconstructed ROOT files following the naming convention:

```text
prad_<RUN>.<FILE>_recon.root
```

For example:

```text
prad_024014.00015_recon.root
```

where:

* `RUN` is a six-digit run number.
* `FILE` is a five-digit file number.

The default input directory is:

```text
/volatile/hallb/prad/wrightso/X17_Random/
```

## Output

Individual analysis outputs are written to:

```text
../rootOutfiles/RandomClusters/<RUN>/
```

Each processed file produces:

```text
TrigVal_<RUN>_<FILE>.root
```

After all jobs complete successfully, the script creates the merged output:

```text
TrigVal_<RUN>_All.root
```

The individual ROOT files are then removed to conserve disk space.

## Requirements

The script requires:

* Python 3
* ROOT
* The ROOT `hadd` utility available in the system path
* The `TrigVal` executable

By default, the executable is expected to be located at:

```text
../TrigVal
```

## Usage

Run the script by specifying the run number:

```bash
python3 Run_triggerStudy1.py <RUN>
```

For example:

```bash
python3 Run_triggerStudy1.py 24014
```

The maximum number of concurrent jobs can be changed using:

```bash
python3 Run_triggerStudy1.py <RUN> --max-parallel <N>
```

For example:

```bash
python3 Run_triggerStudy1.py 24014 --max-parallel 24
```

## Configuration

Several configuration parameters are defined near the top of the script:

| Parameter      | Description                                     |
| -------------- | ----------------------------------------------- |
| `MAX_PARALLEL` | Maximum number of concurrent analysis jobs.     |
| `EXECUTABLE`   | Path to the `TrigVal` executable.               |
| `DATA_DIR`     | Directory containing reconstructed input files. |
| `OUT_BASE`     | Base directory for output ROOT files.           |

These values can be modified to match the local analysis environment.

## Parallel Processing

The script uses Python's `ThreadPoolExecutor` to execute multiple instances of `TrigVal` concurrently. Since each analysis job runs as an independent external process, multiple files can be analyzed simultaneously, substantially reducing the overall processing time for large runs.

## Error Handling

Before processing begins, the script verifies that:

* The input directory exists.
* The `TrigVal` executable is present and executable.
* Input ROOT files matching the requested run are found.

During execution, the script records the exit status and execution time for every analysis job. If any jobs fail, a summary table and the final lines of each job's error output are printed to assist with debugging.

If the `hadd` merge fails, the intermediate ROOT files are preserved to prevent data loss.

## Summary Output

Upon completion, the script reports:

* Run number
* Number of input files processed
* Number of successful and failed jobs
* Total wall-clock execution time
* Location of the merged ROOT file

This summary provides a quick overview of the batch processing results and highlights any files that require further investigation.

## Author Information

This script was originally written by **Rafayel Paremuzyan** ([rafopar@jlab.org](mailto:rafopar@jlab.org)) as `Run_triggerStudy1.py` for automating trigger validation analyses for the Raw Sum Trigger.

It was subsequently adapted and extended by **Erik Wrightson**([wrightso@jlab.org](mailto:wrightso@jlab.org)) for use in this project. The modifications include support for batch processing, configurable parallel execution, automatic merging of per-file ROOT outputs using `hadd`, cleanup of intermediate files, and improved job monitoring and reporting.

---
# 2. Trigger Validation Summary Plotter

## Overview

This Python script reads the output ROOT file produced by the trigger validation analysis and generates a multi-page PDF containing summary plots of trigger performance.

The script recursively searches the ROOT file for histograms, computes trigger efficiencies for multiple trigger thresholds, and produces comparison plots for several trigger quantities. The resulting PDF provides a convenient overview of trigger validation performance for a given dataset.

## Input

The script accepts a single merged ROOT file as input.

For example:

```bash
python3 PlotTriggerValidation.py TrigVal_024014_All.root
```

The input ROOT file is expected to contain the histograms produced by the `TrigVal` analysis.

## Workflow

The script performs the following operations:

1. Opens the specified ROOT file.
2. Recursively searches all directories within the file for histograms.
3. Groups histograms according to their names:

   * Sum
   * Minimum energy
   * Maximum energy
   * XY occupancy
   * Time difference (`deltaT`)
   * Timing correlation
4. Separates each histogram group into:

   * All reconstructed events
   * Events satisfying each VTP trigger threshold
5. Computes trigger efficiencies by dividing the VTP histograms by the corresponding "All" histogram using ROOT's binomial error propagation.
6. Produces a series of comparison plots for each histogram category.
7. Saves all plots into a single multi-page PDF.

## Trigger Efficiency Calculation

For each trigger threshold, the efficiency is calculated as

```text
Efficiency = (Triggered Events) / (All Events)
```

using ROOT's `TH1::Divide()` method with the `"B"` option, which computes binomial uncertainties for each histogram bin.

## Included Plots

The output PDF contains:

* Total energy distributions
* Total energy trigger efficiencies
* Minimum cluster energy distributions
* Minimum cluster energy trigger efficiencies
* Maximum cluster energy distributions
* Maximum cluster energy trigger efficiencies
* Time difference (`ΔT`) distributions
* Time difference trigger efficiencies
* XY occupancy maps for all events and each trigger threshold
* Timing correlation (`Corr`) histograms

For the one-dimensional distributions, each page contains:

* The raw distributions (top panel)
* The corresponding trigger efficiencies (bottom panel)

Several efficiency plots are also provided with a zoomed vertical axis to highlight small efficiency differences near 100%.

## Trigger Thresholds

The script compares the following VTP trigger thresholds:

| Threshold (MeV) |
| --------------: |
|             650 |
|             750 |
|             950 |
|            1100 |
|            1300 |
|            1400 |
|            1500 |

Each threshold is displayed with a unique color throughout the summary plots.

## Output

A multi-page PDF is written in the same directory as the input ROOT file.

The output filename is automatically generated by replacing the ROOT filename with:

```text
<RunNum>AfterMerge.pdf
```

This PDF serves as a compact summary of trigger performance for the analyzed dataset.

## Requirements

The script requires:

* Python 3
* PyROOT
* CERN ROOT

## Usage

Run the script with:

```bash
python3 PlotTriggerValidation.py <MergedROOTFile>
```

Example:

```bash
python3 PlotTriggerValidation.py ../rootOutfiles/RandomClusters/24014/TrigVal_024014_All.root
```

## Notes

* The script recursively searches every directory within the ROOT file, allowing it to process histogram files with nested directory structures.
* Histograms are detached from the input ROOT file after loading so they remain available after the file is closed.
* Histogram grouping is performed using naming conventions (e.g., `Sum`, `Min`, `Max`, `XY`, `deltaT`, and `Corr`), so the input ROOT file should follow the standard `TrigVal` histogram naming scheme.
* The generated PDF is intended as a quick validation report for assessing trigger performance and comparing efficiencies across multiple trigger thresholds.
* Timing correlation (`Corr`) histograms showing the relationship between the trigger time and cluster timing information. These correlations are calculated between the trigger time and either the **first cluster time** or the **last cluster time** in each event, depending on the histogram being displayed.
---
# 3. HyCal XY Trigger Efficiency Analysis

## Overview

This Python script analyzes the spatial dependence of trigger efficiency across the HyCal calorimeter surface. It reads trigger validation output ROOT files, extracts XY occupancy histograms, calculates trigger efficiencies as a function of detector position, and produces a multi-page PDF summarizing the trigger response across HyCal.

The resulting plots allow evaluation of trigger performance versus cluster location, identifying possible geometric effects, acceptance variations, detector inefficiencies, or regions with reduced trigger efficiency.

## Input

The script requires a ROOT file produced by the trigger validation analysis.

Usage:

```bash
python3 PlotXYTriggerEfficiency.py <ROOT_FILE>
```

Example:

```bash
python3 PlotXYTriggerEfficiency.py ../rootOutfiles/RandomClusters/24014/TrigVal_024014_All.root
```

The input file is expected to contain XY histograms following the naming convention produced by `TrigVal`.

## Workflow

The script performs the following steps:

1. Opens the provided ROOT file.
2. Recursively searches through all directories to find histograms.
3. Identifies XY position histograms for:

   * All events
   * Each VTP trigger threshold
4. Removes the central HyCal region from efficiency calculations.
5. Calculates trigger efficiencies as a function of detector position.
6. Produces two-dimensional efficiency maps and one-dimensional projections.
7. Saves all plots into a multi-page PDF.

## Histogram Selection

The script searches for histograms containing:

```text
XY
```

These histograms are separated into:

### All Events

The reference distribution containing all reconstructed events:

```text
XY_All
```

### Triggered Events

The distributions for each VTP trigger threshold:

```text
XY_VTP
```

The trigger thresholds analyzed are:

| VTP Threshold (MeV) |
| ------------------: |
|                 650 |
|                 750 |
|                 950 |
|                1100 |
|                1300 |
|                1400 |
|                1500 |

## Trigger Efficiency Calculation

The trigger efficiency is calculated independently for each HyCal position bin:

```text
Efficiency(x,y) = Triggered(x,y) / All(x,y)
```

The calculation uses ROOT's binomial uncertainty option:

```cpp
TH1::Divide(..., "B")
```

which provides statistical uncertainties appropriate for efficiency measurements.

## Central Region Removal

The central HyCal region is excluded from the efficiency calculation because it corresponds to the beam hole or otherwise unusable detector region.

The excluded bins are:

```text
X bins: 16–19
Y bins: 16–19
```

These bins are removed from both the numerator and denominator before calculating projected efficiencies.

## Output

The script creates a PDF file:

```text
_XY_Ratios.pdf
```

in the same directory as the input ROOT file.

The PDF contains:

### XY Occupancy Maps

Two-dimensional detector maps showing:

* All reconstructed events
* Events passing each VTP trigger threshold

The maps are displayed using logarithmic color scaling to improve visualization of regions with different event populations.

### XY Trigger Efficiency Maps

For each trigger threshold, the script produces:

* A two-dimensional efficiency map:

```text
Efficiency(x,y)
```

showing trigger efficiency across the HyCal face.

* An X projection:

```text
Efficiency(x)
```

* A Y projection:

```text
Efficiency(y)
```

These projections provide a clearer view of position-dependent trigger behavior along each detector axis.

## Interpreting the Results

The XY efficiency maps can be used to investigate:

* Uniformity of trigger response across HyCal.
* Inefficient detector regions.
* Effects near calorimeter boundaries.
* Trigger threshold dependence.
* Possible geometric or electronics-related inefficiencies.

A uniform efficiency map near 100% indicates consistent trigger behavior, while localized decreases may indicate detector acceptance effects or trigger logic limitations.

## Requirements

The script requires:

* Python 3
* CERN ROOT
* PyROOT

## Notes

* The script assumes the input ROOT file follows the histogram naming conventions generated by the `TrigVal` trigger validation analysis.
* Histograms are detached from the ROOT file before closing so they remain available for plotting.
* Efficiency calculations use statistical binomial errors.
* The script focuses on spatial trigger performance and complements timing and energy-based trigger validation studies.

## Relationship to Trigger Validation Analysis

This script is designed to be run after the trigger validation processing workflow. The merged ROOT output from the batch trigger validation analysis provides the required XY histograms used to evaluate the trigger efficiency across the HyCal detector.
