# Flask server
#   https://flask.palletsprojects.com/en/2.2.x/
#   python3 flask_api.py

import json

from flask import Flask, request
from flask_cors import CORS
from functions import (
    __generate_spectrum,
    __param_check,
    __process_background,
    __process_spectrum,
    __find_peaks,
)

app = Flask(__name__)
CORS(app)


@app.route("/", methods=["GET"])
def ftir():
    return "<h1 style='color:blue'>Raston Lab FTIR API</h1>"


@app.route("/spectrum", methods=["POST"])
def spectrum():
    # put incoming JSON into a dictionary
    params = json.loads(request.data)

    # verify user input is valid
    if not __param_check(params):
        return {
            "success": False,
            "text": "Parameter check failed",
        }

    # perform:
    #   --> transmission spectrum of gas sample (calc_spectrum)
    spectrum, error, message = __generate_spectrum(params)
    if error:
        return {
            "success": False,
            "text": message,
        }

    # perform:
    #   --> blackbody spectrum of source (sPlanck)
    #   --> transmission spectrum of beamsplitter and cell windows
    #   --> detector response spectrum
    processed_spectrum = __process_spectrum(params, spectrum)

    # https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.get
    x_value, y_value = processed_spectrum.get("transmittance_noslit")

    # convert dictionary values to strings and return as JSON
    return {
        "success": True,
        "x": list(x_value),
        "y": list(map(str, y_value)),
    }


@app.route("/background", methods=["POST"])
def background():
    try:
        # put incoming JSON into a dictionary
        data = json.loads(request.data)
    
        # verify user input is valid
        if not __param_check(data):
            return {
                "success": False,
                "text": "Parameter check failed",
            }
    
        # perform:
        #   --> transmission spectrum of gas sample (calc_spectrum)
        spectrum, error, message = __generate_spectrum(data)
        if error:
            return {
                "success": False,
                "text": message,
            }
    
        # perform:
        #   --> set all y-values to one
        try:
            background_spectrum = __process_background(spectrum)
        except:
            return {
                "success": False,
                "text": "Background Failure"
            }
    
        # perform:
        #   --> blackbody spectrum of source (sPlanck)
        #   --> transmission spectrum of beamsplitter and cell windows
        #   --> detector response spectrum
        processed_spectrum = __process_spectrum(data, background_spectrum, True)
    
        if processed_spectrum is None:
            return {
                "success": False,
                "text": "Issue Processing Data"
            }
    
        # https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.get
        x_value, y_value = processed_spectrum.get("transmittance_noslit")
    
        # convert dictionary values to strings and return as JSON
        return {
            "success": True,
            "x": list(x_value),
            "y": list(map(str, y_value)),
        }

    # perform:
    #   --> blackbody spectrum of source (sPlanck)
    #   --> transmission spectrum of beamsplitter and cell windows
    #   --> detector response spectrum
    processed_spectrum = __process_spectrum(data, background_spectrum)

    if processed_spectrum is None:
        return {
            "success": False,
            "text": "Issue Processing Data"
        }
      
    except:
       return {
            "success": False,
            "text": "Issue Processing Data"
        } 


@app.route("/find_peaks", methods=["POST"])
def find_peaks():
    data = json.loads(request.data)

    peaks = __find_peaks(
        data["x"],
        data["y"],
        float(data["lowerbound"]),
        float(data["upperbound"]),
        float(data["threshold"]),
    )

    if peaks:
        return {"success": True, "peaks": peaks, "text": None}
    else:
        return {"success": False, "peaks": None, "text": "Failed to find peaks"}


# set debug to false in production environment
if __name__ == "__main__":
    app.run(host="0.0.0.0")
