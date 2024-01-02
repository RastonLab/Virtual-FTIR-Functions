# Virtual-FTIR-Functions
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/RastonLab/Virtual-FTIR-Functions/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/RastonLab/Virtual-FTIR-Functions/tree/main)

This repository contains the back-end of the [Virtual FTIR Spectrometer](https://github.com/RastonLab/Virtual-FTIR-Spectrometer).

## Flask

- `wsgi.py`

  - This file is used for deployment by Gunicorn and Nginx on the server.

- `app.py`

  - Contains a program that accepts API calls for both spectrum generation and background sample generation. This program accepts JSON through a `fetch()` or `cURL` command and returns JSON of X and Y coordinates.

- `processing.py`

  - This module contains four functions: `process_spectrum`, `generate_background`, `generate_spectrum`, and `find_peaks`.  These functions are used for the generation and anaylsis of realistic spectra. Find more details about each function [here](#processingpy-functions).

- `processing_utils.py`

  - This utility module contains helper functions for the `processing.py` module. These functions include functions to approxiamate physical components of an FTIR Spectrometer, approximate realistic noise, check user parameters, and calculate the resolution of the spectra. Find more details about the individual functions [here](#processing_utilspy-functions).

## Scripts

- `download_hitran.py`

  - This script downloads HITRAN data files locally. This is useful so the Flask application doesn't need to download data files during user queries.

## Installation

Information on how to run the back-end can be found in the repository's [wiki page](../../wiki).

> **NOTE**: Make sure in the frontend that `src/dictionaries/constants.js` has the proper fetch URLs uncommented.

## The Process - Breakdown of App.py

**A note on vocabulary:** in this documentation "_ideal spectrum_" means either a background or sample spectrum with no noise and is unaffected by the physical components of the spectrometer.

<!-- TODO describe the parameters for each POST Request -->
There are two sets of parameters the server can receive based on what request is being made. If the request is for Spectrum or Background, the parameters will include:

```
  beamsplitter
  detector
  window
  source
  resolution
  zeroFill
  scan
  mole
  molecule
  pressure
  waveMax
  waveMin
```
`beamsplitter`, `detector`, `window`, and `source` are used to select which [component spectra](#component-functions) will be used in `process_spectra` ([details](#process_spectrum)). `resolution` and `zeroFill` are used by `calc_wstep` to determine the resolution of the spectrum ([details](#calc_wstep)). `scan` is used by `multiscan` and determines the number of scans to be simulated ([details](#multiscan)). The rest of the parameters (`mole`, `molecule`, `pressure`, `waveMax`, `waveMin`) are used in `generate_spectrum` to create the spectrum ([details](#generate_spectrum)). 

For Find Peaks requests, the parameters will include:

```
  x-data
  y-data
  lowerbound
  upperbound
  threshold
```
`x-data` and `y-data` contain the x and y-values for the spectrum to analyze. `lowerbound` gives the lowest x-value in the range to analyze; `upperbpund` gives the higest x-value in the range to analyze. `threshold` gives the lowest y-value to consider a 'peak' in the data.

- Sample

  - Sample recives a list of parameters from the POST Request in the form of a JSON. After running some basic [parameter checking](#param_check), Sample calls `generate_spectrum`. This function returns an ideal spectrum in the form of a RADIS Spectrum object based on the user selected parameters. This ideal spectrum is then passed to the `process_spectrum` function which alters the spectrum to more closly resemble a realistic spectrum. For more detail on this function, go [here](#process_spectrum). Then the x and y values from the resulting spectrum are sent back to the user.

- Background

  - Background takes in parameters, checks them, and calls `generate_spectrum` the same as Spectrum. After these preliminary steps, Background then calls `generate_background` to get a spectrum where all the y-values are 1 ([details](#generate_background)). By calling `generate_spectrum` and then calling `generate_background` with the resulting spectrum, the x-values and overall resolution of the spectrum is accurate; the background data lines up with a sample spectrum with the same parameters. From here, the spectrum is sent to `process_spectrum` so it will more closely resemble a real spectrum. Like Spectrum, this function sends the x and y-values of the resulting spectrum back to the user.

- Find Peaks

  - Find Peaks simply calls the `find_peaks` function using the parameters  the user selected. For more details on `find_peaks`, see the [`find_peaks`](#find_peaks) breakdown.

## Function Details

### Processing.py Functions

#### `process_spectrum`

  - This function takes an ideal spectrum and returns a spectrum that is approximately the spectrum that would be generated by a physical spectrometer. This is achieved by creating additional spectra based on mathematical functions that approximate the behavior of FTIR components. The user chooses which of these spectra will be used by picking a [Beamsplitter](#beamsplitters), [Cell Window](#cell-windows), [Detector](#detector), and a source (Globar or Tungsten) for the [Blackbody Spectrum](#blackbody). Those spectra are then multiplied into the base spectrum and realistic noise is added. Find more detail on the component functions [here](#component-functions) and find more detail about noise generation [here](#multiscan).

#### `generate_background`

  - This function takes the wavenumbers (x-values) of a spectra and sets all y-values to one to simulate an ideal background spectrum.

#### `generate_spectrum`

  - This function utilizes [RADIS](https://radis.readthedocs.io/en/latest/index.html) and the `calc_spectrum` function provided ([details](https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum)) to calculate an ideal spectrum based on the [parameters provided](#the-process---breakdown-of-apppy) in the `params` parameter and the resolution determined by `calc_wstep` ([details](#calc_wstep)).

#### `find_peaks`

  - This function takes the x and y-values of a spectrum and finds all the peaks in the data using [RADIS tools](https://specutils.readthedocs.io/en/stable/api/specutils.fitting.find_lines_threshold.html#specutils.fitting.find_lines_threshold). Then the desired (emmission) peaks are returned from that data based on the minimum value for a y-value to be considered a peak (`threshold`).

---

### Processing_utils.py Functions

#### `zeroY`

  - This function sets the y-values for the given spectrum to 1 in order to emulate a background sample.

#### `param_check`

  - This function serves as a sanity check for the parameters sent in by the user. They are thoroughly checked on the frontend before the post request to this server is made.

#### `calc_wstep`

  - This function calculates the appropriate spectrum resolution/wstep based on the user given parameters `resolution` and `zero_fill`. This function is based off of scientific data collected by the RastonLab.

#### `multiscan`

  - The purpose of this function is to add noise to the given spectra while limiting the amount of memory space used at any given time. The amount of noise added in total is determined by the `num_scans` variable. The amount of noise added in a particular iteration is determined by the value of `scans_per_group`, which is currently set to `10`. In each iteration, a 2d array is created containing random numbers determined by `np.random.normal`; in this 2d array there is a column for each x-value in the provided spectra and a row for each scan being simulated. All of these values are added together and then divided by the total number of scans (for normalization purposes). The resulting 1d array is then added into the spectrum.

#### `get_component_spectra`

  - This function is a helper function that generates all the the component spectra. It was created to reduce the amount of code in the `process_spectrum`. Additionally, this function normalizes the blackbody spectrum and both detector spectra. For more detail on this, see [Component Functions](#component-functions).

#### Component Functions

 These functions were designed to emulate several physical components within a real FTIR Spectrometer and their effect on spectra. This was accomplished by analyzing the data collected by the RastonLab depicting how each physical component alters the spectra. Then a mathematical function was created for each component through [curve fitting](https://en.wikipedia.org/wiki/Curve_fitting#:~:text=Curve%20fitting%20is%20the%20process,points%2C%20possibly%20subject%20to%20constraints.). These functions are then used to create y-values that correspond to the list of x-values from the ideal spectrum. These x-y value pairs are then used to make Spectrum objects for each component (`get_component_spectra`).

 Following are the graphs associated with each function. Some functions may have more details as necessary.

 **NOTE:** In many of these functions, there is a variable named `x_um` that is set equal to `10000 / spectrum`. This (and any other instance of `10000 / spectrum`) is there to handle unit conversions.

#### Blackbody

  - #### `__sPlanck`

    This function takes a parameter called `source_temp`. This temperature is either `1200`, which is associated with Globar, or `3400`, which is associated with Tungsten. The spectrum that is generated with this function is normalized to 1 using the RADIS [normalize function](https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.normalize)

    - Globar
    
      ![Blackbody-Globar](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/689da3ff-5637-484f-9cd8-6bb5497e86ad)

    - Tungsten
   
      ![Blackbody-Tungsten](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/bda16a62-1397-458a-8eb6-8f655e647c0e)

#### Cell Windows

  - #### `__CaF2`

    ![Window-CaF2](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/a8a18cbb-a31a-4432-acaa-78451c4c5cfd)

  - #### `__ZnSe`
 
    ![Window-ZnSe](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/3b9147eb-a6cc-4342-a26b-1afec82bdff2)

#### Beamsplitters

  - #### `AR_ZnSe`
 
    ![Beamsplitter-AR_ZnSe](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/1570ba60-44a1-4afe-92e1-a003704bb256)

  - #### `AR_CaF2`
 
    ![Beamsplitter-AR_CaF2](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/839c9b70-7c14-42f4-81f5-d4a0a98bf205)

#### Detector

  - #### `__InSb`

    The spectrum generated by this function is normalized to 2 using the RADIS [normalize function](https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.normalize), multiplying the resulting y-values by 2, and then making a new spectrum using those new y-values.

    When this Detector is chosen, the spectrum for `__sapphire` is also multiplied into the spectrum.
 
    ![Detector-InSb](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/461d479f-e6be-4e49-841a-60f630a70c5b)

    - #### `__sapphire`

      **NOTE:** This function takes a large number to the `6903.57039` power. The size of this exponent is unavoidable due to the curve fitting technique used to generate this function. As a result, this function invokes an `Overflow Error` at runtime. There have not been any noticeable problems caused by this error.
   
      ![Sapphire](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/673cbf31-2bd7-4949-9954-7dcca95881f2)

  - #### `__MCT`

    The spectrum generated with this function is normalized to 1 using the RADIS [normalize function](https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.normalize)

    When this Detector is chosen, the spectrum for `__ZnSe` is also mulitiplied into the spectrum, regardless of whether or not `__ZnSe` was multiplied in earlier.
 
    ![Detector-MCT](https://github.com/Brennaser/Virtual-FTIR-Functions/assets/54820278/7ad8a7ed-90cd-4f92-a869-79b8265f7480)
