# Virtual-FTIR-Functions

This repository contains the back-end of the [Virtual FTIR Spectrometer](https://github.com/RastonLab/Virtual-FTIR-Spectrometer).

## [api](https://github.com/RastonLab/Virtual-FTIR-Functions/tree/main/api)

This directory contains

- `background.py`
  - Contains a program that strictly generates a background sample. This program uses `data.json` for its input and is capable of graphing using `mathplotlib`.
- ## `fast_api.py`
  - Contains a program that accepts API calls for both a background sample and spectrum. This program accepts JSON through a `fetch()` or `cURL` command and returns JSON of X and Y coordinates.
- ## `flask_api`
  - Contains a program that accepts API calls for both a background sample and spectrum. This program accepts JSON through a `fetch()` or `cURL` command and returns JSON of X and Y coordinates.
- ## `functions.py`
  - Contains functions used in all other files. Designed so edits only need to be made to one file for main functionality updates and changes.
- ## `spectrum.py`
  - Contains a program that strictly generates a spectrum. This program uses `data.json` for its input and is capable of graphing using `mathplotlib`.

## [jupyter-notebook](https://github.com/RastonLab/Virtual-FTIR-Functions/tree/main/jupyter-notebook)

This directory holds a Jupyter Notebook that was used in the testing and development of the spectrum modification functions. These functions emulate the blackbody spectrum (sPlanck), the beamsplitter, the cell windows, and the detector response spectrum.

## [create-virtual-environments.sh](https://github.com/RastonLab/Virtual-FTIR-Functions/blob/main/create-virtual-environments.sh)

This script creates a virtual python environments (venv) needed to run the files in the `api` and `jupyter-notebook` directories. The script uses `requirements.txt` to install the required dependencies to a venv that is used for both directories.

To activate the virtual environment (venv), run `source venv/bin/activate` in the terminal. To deactivate the virtual environment (venv), run `deactivate` in the terminal.
