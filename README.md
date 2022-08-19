# Virtual-FTIR-Functions

This repository contains the back-end of the [Virtual FTIR Spectrometer](https://github.com/RastonLab/Virtual-FTIR-Spectrometer).

## [flask](https://github.com/RastonLab/Virtual-FTIR-Functions/tree/main/flask)

This directory holds the Flask web framework that the front-end reaches out to.

## [jupyter-notebook](https://github.com/RastonLab/Virtual-FTIR-Functions/tree/main/jupyter-notebook)

This directory holds a Jupyter Notebook that was used in the testing and development of the spectrum modification functions. These functions emulate the blackbody spectrum (sPlanck), the beamsplitter, the cell windows, and the detector response spectrum.

## [create-virtual-environments.sh](https://github.com/RastonLab/Virtual-FTIR-Functions/blob/main/create-virtual-environments.sh)

This script creates a virtual python environments (venv) needed to run the files in the `flask` and `jupyter-notebook` directories. The script uses `requirements.txt` to install the required dependencies to a venv that is used for both directories.

To activate the virtual environment (venv), run `source venv/bin/activate` in the terminal. To deactivate the virtual enviornment (venv), run `deactivate` in the terminal.
