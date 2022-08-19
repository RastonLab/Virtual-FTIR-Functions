from flask import Flask, request
from flask_cors import CORS

import json

from functions import __param_check, __generate_spectra, __generate_background

app = Flask(__name__)
CORS(app)


@app.route("/", methods=["GET"])
def ftir():
    return "<h1>ftir</h1>"


@app.route("/fetch_background", methods=["POST"])
def fetch_background():
    # put incoming JSON into a dictionary
    data = json.loads(request.data)

    print(data)

    # verify the information in the dictionary
    __param_check(data)

    print("----- output verified params to console as self-check -----")
    for key, value in data.items():
        print("  %s: %s" % (key, value))

    # perform: transmission spectrum of gas sample (calc_spectrum)
    #      --> blackbody spectrum of source (sPlanck)
    #      --> transmission spectrum of beamsplitter and cell windows
    #      --> detector response spectrum
    print("----- start __generate_background() -----")
    result = __generate_background(data)
    print("----- end __generate_background() -----")

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

@app.route("/post_json", methods=["POST"])
def process_json():
    # put incoming JSON into a dictionary
    data = json.loads(request.data)

    print(data)

    # verify the information in the dictionary
    __param_check(data)

    print("----- output verified params to console as self-check -----")
    for key, value in data.items():
        print("  %s: %s" % (key, value))

    # perform: transmission spectrum of gas sample (calc_spectrum)
    #      --> blackbody spectrum of source (sPlanck)
    #      --> transmission spectrum of beamsplitter and cell windows
    #      --> detector response spectrum
    print("----- start __generate_spectra() -----")
    result = __generate_spectra(data)
    print("----- end __generate_spectra() -----")

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
