# Flask server
#   https://flask.palletsprojects.com/en/2.2.x/
#   python3 flask_api.py

import json

from flask import Flask, request
from flask_cors import CORS
from processing import (
    generate_spectrum,
    generate_background,
    process_spectrum,
    find_peaks
)
from processing_utils import param_check

app = Flask(__name__)
CORS(app)
try:
	with open("version.txt","r") as f:
		version = f.read()
		app.config["VERSION"] = version
except:
     print("no version file found")

@app.route("/", methods=["GET"])
def ftir() -> str:
    if "VERSION" not in app.config:
      app.config["VERSION"] = "0.0.0"
    return "<h1 style='color:blue'>Raston Lab FTIR API%s</h1>" % (" - Version "+app.config["VERSION"])

@app.route("/sample", methods=["POST"])
def sample() -> dict[bool, list[float], list[float]]:
    # put incoming JSON into a dictionary
    params = json.loads(request.data)

    # verify user input is valid
    if not param_check(params):
        return {
            "success": False,
            "text": "One of the given parameters was invalid. Please change some settings and try again.",
        }

    # perform:
    #   --> transmission spectrum of gas sample (calc_spectrum)
    spectrum, error, message = generate_spectrum(params)
    if error:
        return {
            "success": False,
            "text": message,
        }

    # perform:
    #   --> blackbody spectrum of source (sPlanck)
    #   --> transmission spectrum of beamsplitter and cell windows
    #   --> detector response spectrum
    processed_spectrum = process_spectrum(params, spectrum)

    # https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.get
    x_value, y_value = processed_spectrum.get("transmittance_noslit")

    # convert dictionary values to strings and return as JSON
    return {
        "success": True,
        "x": list(x_value),
        "y": list(map(str, y_value)),
    }


@app.route("/background", methods=["POST"])
def background() -> dict[bool, list[float], list[float]]:
    # put incoming JSON into a dictionary
    data = json.loads(request.data)

    # verify user input is valid
    if not param_check(data):
        return {
            "success": False,
            "text": "One of the given parameters was invalid. Please change some settings and try again.",
        }

    # perform:
    #   --> transmission spectrum of gas sample (calc_spectrum)
    spectrum, error, message = generate_spectrum(data)
    if error:
        return {
            "success": False,
            "text": message,
        }

    # perform:
    #   --> set all y-values to one
    try:
        background_spectrum = generate_background(spectrum)
    except:
        return {
            "success": False,
            "text": "There was an issue collecting the background spectra."
        }

    # perform:
    #   --> blackbody spectrum of source (sPlanck)
    #   --> transmission spectrum of beamsplitter and cell windows
    #   --> detector response spectrum
    processed_spectrum = process_spectrum(data, background_spectrum)
    
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


@app.route("/find_peaks", methods=["POST"])
def handle_peaks() -> dict[bool, dict[float, float], str]:
    data = json.loads(request.data)

    peaks, error = find_peaks(
        data["x"],
        data["y"],
        float(data["threshold"]),
    )

    if (peaks): 
        return {"success": True, "peaks": peaks, "text": error}

    return {"success": False, "peaks": peaks, "text": error}


# set debug to false in production environment
if __name__ == "__main__":
    app.run(host="0.0.0.0")
