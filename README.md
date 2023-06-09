# Virtual-FTIR-Functions
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/RastonLab/Virtual-FTIR-Functions/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/RastonLab/Virtual-FTIR-Functions/tree/main)

This repository contains the back-end of the [Virtual FTIR Spectrometer](https://github.com/RastonLab/Virtual-FTIR-Spectrometer).

## Flask

- `app.py`

  - Contains a program that accepts API calls for both spectrum generation and background sample generation. This program accepts JSON through a `fetch()` or `cURL` command and returns JSON of X and Y coordinates.

- `functions.py`

  - Contains functions used in all other files. Designed so edits only need to be made to one file for main functionality updates and changes.

- `wsgi.py`

  - This file is used for deployment by Gunicorn and Nginx on the server.

## scripts

- `download_hitran.py`

  - This script downloads HITRAN data files locally. This is useful so the Flask application doesn't need to download data files during user queries.

- `virtual_environment.sh`

  - This script creates a virtual python environments (venv) needed to run `app.py`. The script uses `requirements.txt` to install the required dependencies to a venv that is used for both directories. To activate the virtual environment (if your editor or terminal does not automatically activate it), run `source venv/bin/activate`. To deactivate the virtual environment, run `deactivate` in the terminal.
    - NOTE: This script is designed to be run from outside the `scripts` directory.

---

## How to run

1. Create and setup the Python virtual environment:

   ```
   ./script/virtual_environment.sh
   ```

2. Activate the virtual environment:

   ```
   source venv/bin/activate
   ```

3. To run the Flask application:

   ```
   python app.py
   ```

   or

   ```
   python wsgi.py
   ```

**NOTE**: Make sure in the frontend that `src/routes/ExperimentalSetup.jsx` has the proper fetch URLs uncommented: https://github.com/RastonLab/Virtual-FTIR-Spectrometer/blob/main/src/routes/ExperimentalSetup.jsx

---

## Flask Test Query

- `cURL`

  - This cURL command has be used in a terminal to test the applications ability to properly return the X and Y coordinates.

  ```
  curl -X POST localhost:5000/spectrum \
      -H "Content-type: application/json" \
      -d "{ \
          \"minWave\" : 1900, \
          \"maxWave\" : 2300, \
          \"molecule\" : \"CO\", \
          \"pressure\" : 0.01, \
          \"resolution\" : 1, \
          \"numScan\" : 1, \
          \"zeroFill\" : 0, \
          \"source\" : 3100, \
          \"beamsplitter\" : \"AR_ZnSe\", \
          \"cellWindow\" : \"CaF2\", \
          \"detector\" : \"MCT\" \
      }"
  ```
