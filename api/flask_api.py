# Flask server
#   https://flask.palletsprojects.com/en/2.2.x/
#   python3 flask_api.py

import json

from flask import Flask, request
from flask_cors import CORS
from functions import __generate_background, __generate_spectra, __param_check

app = Flask(__name__)
CORS(app)


@app.route("/", methods=["GET"])
def ftir():
    return "<h1 style='color:blue'>Raston Lab FTIR API</h1>"


@app.route("/background", methods=["POST"])
def fetch_background():
    # put incoming JSON into a dictionary
    data = json.loads(request.data)

    # verify user input is valid
    if not __param_check(data):
        return {
            "success": False,
            "text": "Parameter check failed",
        }

    # perform: transmission spectrum of gas sample (calc_spectrum)
    #      --> blackbody spectrum of source (sPlanck)
    #      --> transmission spectrum of beamsplitter and cell windows
    #      --> detector response spectrum
    result = __generate_background(data)

    if result != False:
        # convert dictionary values to strings and return as JSON
        return {
            "success": True,
            "x": list(result.keys()),
            "y": [str(flt) for flt in result.values()],
        }
    else:
        return {
            "success": False,
            "text": "No line in the specified wavenumber range",
        }


@app.route("/spectra", methods=["POST"])
def process_json():
    # put incoming JSON into a dictionary
    data = json.loads(request.data)

    # verify user input is valid
    if not __param_check(data):
        return {
            "success": False,
            "text": "Parameter check failed",
        }

    # perform: transmission spectrum of gas sample (calc_spectrum)
    #      --> blackbody spectrum of source (sPlanck)
    #      --> transmission spectrum of beamsplitter and cell windows
    #      --> detector response spectrum
    result = __generate_spectra(data)

    if result != False:
        # convert dictionary values to strings and return as JSON
        return {
            "success": True,
            "x": list(result.keys()),
            "y": [str(flt) for flt in result.values()],
        }
    else:
        return {
            "success": False,
            "text": "No line in the specified wavenumber range",
        }


if __name__ == "__main__":
    app.run(debug=True)
