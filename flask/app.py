from flask import Flask, request
from flask_cors import CORS
from radis import calc_spectrum

import json
import math
import numpy as np

app = Flask(__name__)
CORS(app)


@app.route("/", methods=["GET"])
def ftir():
    return "<h1>ftir</h1>"

@app.route("/fetch_background", methods=["POST"])
def fetch_background():
    data = json.loads(request.data)
    # Need to add Param Check
    print(calc_spectrum(
            data["minWave"],
            data["maxWave"],
            molecule= None,
            # isotope="1,2,3",
            # pressure=data["pressure"],
            Tgas=294.15,  # hardcode
            path_length=10,  # hardcode
            wstep=0.5,  # (cm^-1)
            verbose=False,  # hides HITRAN output
            databank="hitran",
            warnings={"AccuracyError": "ignore"},
        ))
    result = ""

    # convert dictionary values to strings and return as JSON
    return {
        "success": True,
        "x": list(result.keys()),
        "y": [str(flt) for flt in result.values()],
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


# --------------------------------------
# ---- spectra calculation functions ----
# --------------------------------------


def __KBr(data):
    if data == None:
        return False

    for x in data:
        datapoint = (0.92267) / (1 + (25.66477 / (x / 1000)) ** -12.35159) ** 0.17344
        data[x] = datapoint * data[x]

    return data


def __CaF2(data):
    if data == None:
        return False

    for x in data:
        datapoint = (0.93091) / (1 + (11.12929 / (x / 1000)) ** -12.43933) ** 4.32574
        data[x] = datapoint * data[x]

    return data


def __ZnSe(data):
    if data == None:
        return False

    for x in data:
        x_um = x / 1000
        datapoint = (0.71015) / (
            (1 + (20.99353 / x_um) ** -19.31355) ** 1.44348
        ) + -0.13265 / (2.25051 * math.sqrt(math.pi / (4 * math.log(2)))) * math.exp(
            -4 * math.log(2) * ((x_um - 16.75) ** 2) / (2.25051**2)
        )
        data[x] = datapoint * data[x]

    return data


def __sapphire(data):
    if data == None:
        return False

    for x in data:
        # Gets accurate graph with numpy float128 but throws runtime overflow error
        # datapoint = Decimal(0.78928) / Decimal(1 + (11.9544 / (x / 1000)) ** -12.07226 ) ** (Decimal(6903.57039))
        datapoint = np.float128(0.78928) / np.float128(
            1 + (11.9544 / (x / 1000)) ** -12.07226
        ) ** (np.float128(6903.57039))
        dp2 = datapoint * np.float128(data[x])
        data[x] = dp2

    return data


def __AR_ZnSe(data):
    if data == None:
        return False

    for x in data:
        x_um = x / 1000
        datapoint = (
            (0.82609) / ((1 + ((34.63971 / x_um) ** -8.56269)) ** 186.34792)
            + -0.47
            / (0.55 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 1.47) ** 2) / (0.55**2))
            + -0.03456
            / (0.4 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 2.88) ** 2) / (0.4**2))
            + -0.009
            / (0.3 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 6.16) ** 2) / (0.3**2))
            + -0.09
            / (1 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 16.2) ** 2) / (1**2))
            + -0.08
            / (1 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 17.4) ** 2) / (1**2))
            + 1.12
            / (8 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 9.5) ** 2) / (8**2))
            + 0.11546
            / (2 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 4.9) ** 2) / (2**2))
            + 0.21751
            / (2 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 2.6) ** 2) / (2**2))
            + -0.05
            / (0.07 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 0.8) ** 2) / (0.07**2))
        )

        data[x] = datapoint * data[x]

    return data


def __AR_CaF2(data):
    if data == None:
        return False

    for x in data:
        x_um = x / 1000
        datapoint = (
            (0.9795) / ((1 + ((18.77617 / x_um) ** -6.94246)) ** 91.98745)
            + -0.06
            / (0.08 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 0.76) ** 2) / (0.08**2))
            + -0.06
            / (0.2 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * (x_um - 1.06) ** 2 / 0.20**2)
            + -0.6
            / (3.0 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 4.85) ** 2) / (3.0**2))
            + -0.35
            / (1.0 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 9.40) ** 2) / (1.00**2))
            + 0.05
            / (0.8 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 2.60) ** 2) / (0.8**2))
            + 0.04
            / (0.5 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 7.75) ** 2) / (0.50**2))
            + -0.01
            / (0.6 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 6.55) ** 2) / (0.6**2))
            + 0.01
            / (0.5 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 1.82) ** 2) / (0.5**2))
        )
        data[x] = datapoint * data[x]

    return data


def __InSb(data):
    if data == None:
        return False

    for x in data:
        x_um = x / 1000
        datapoint = 1.97163e11 * (1 / (1 + math.exp(-(x_um - 5.3939) / 1.6624))) * (
            1 - 1 / (1 + math.exp(-(x_um - 5.3939) / 0.11925))
        ) + (3.3e10) / (2.44977 * math.sqrt(math.pi / (4 * math.log(2)))) * math.exp(
            -4 * math.log(2) * ((x_um - 5) ** 2) / (2.44977**2)
        )
        data[x] = (datapoint + np.random.normal(0, 200000000)) * data[x]

    return data


def __MCT(data):
    if data == None:
        return False

    for x in data:
        x_um = x / 1000
        datapoint = (
            (1.98748 * (10**9))
            + (2.10252 * (10**10))
            * (1 / (1 + math.exp(-(x_um - 20.15819) / 5.73688)))
            * (1 - 1 / (1 + math.exp(-(x_um - 20.15819) / 1.11659)))
            + (1.3 * (10**9))
            / (2 * math.sqrt(math.pi / (4 * math.log(2))))
            * math.exp(-4 * math.log(2) * ((x_um - 18.6) ** 2) / (2**2))
        )
        data[x] = (datapoint + np.random.normal(0, 20000000)) * data[x]

    return data


def __sPlanck(spectrum, temp):
    H = 6.62606957e-34
    C = 2.99792458e8
    K_B = 1.3806488e-23

    if spectrum == None:
        return False

    for x in spectrum:
        x2 = x * (10**-9)
        p = ((0.2 * H * (C**2)) / ((x2**4) * x)) * (
            1 / (math.exp((H * C) / (x2 * K_B * temp)) - 1)
        )
        spectrum[x] = spectrum[x] * p

    return spectrum


# -------------------------------------
# ---------- helper functions ----------
# ------------------------------------
def __loadData(s):
    data = {}

    for key, val in zip(s[0], s[1]):
        data[float(key)] = float(val)

    return data


def __param_check(data):
    print("----- check number of keys -----")
    if len(data) == 11:
        print("  good!")
    else:
        print("  not enough params. total params: %s" % (len(data)))
        return "not enough params. total params: %s" % (len(data))

    print("----- check if params are correct -----")
    valid_params = [
        "minWave",
        "maxWave",
        "molecule",
        "pressure",
        "resolution",
        "numScan",
        "zeroFill",
        "source",
        "beamsplitter",
        "cellWindow",
        "detector",
    ]
    for key, value in data.items():
        if key in valid_params:
            if (data[key] == "") or (data[key] == None):
                print("  error with key: %s. Value is: %s" % (key, value))
                return "Error with key: %s. Value is: %s" % (key, value)
            else:
                print("  %s: %s" % (key, value))
        else:
            print("  error with key: %s. Value is: %s" % (key, value))
            return "Error with key: %s. Value is: %s" % (key, value)


def __calc_wstep(resolution, zero_fill):
    
    wstep = 0

    if resolution == 1:

        if zero_fill == 0:
            wstep = 0.481927711
        elif zero_fill == 1:
            wstep = 0.240963855
        elif zero_fill == 2:
            wstep = 0.120481928

    elif resolution == 0.5:

        if zero_fill == 0:
            wstep = 0.240963855
        elif zero_fill == 1:
            wstep = 0.120481928
        elif zero_fill == 2:
            wstep = 0.060240964

    elif resolution == 0.25:

        if zero_fill == 0:
            wstep = 0.120481928
        elif zero_fill == 1:
            wstep = 0.060240964
        elif zero_fill == 2:
            wstep = 0.030120482

    elif resolution == 0.125:

        if zero_fill == 0:
            wstep = 0.060240964
        elif zero_fill == 1:
            wstep = 0.030120482
        elif zero_fill == 2:
            wstep = 0.015060241

    elif resolution == 0.0625:

        if zero_fill == 0:
            wstep = 0.030120482
        elif zero_fill == 1:
            wstep = 0.015060241
        elif zero_fill == 2:
            wstep = 0.00753012

    return wstep

def __generate_spectra(data):
    try:
        wstep = __calc_wstep(data["resolution"], data["zeroFill"])
        # ----- a.) transmission spectrum of gas sample -----
        # https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum
        s = calc_spectrum(
            data["minWave"],
            data["maxWave"],
            molecule=data["molecule"],
            isotope="1,2,3",
            pressure=data["pressure"],
            Tgas=294.15,  # hardcode
            path_length=10,  # hardcode
            wstep=wstep,  # (cm^-1)
            verbose=False,  # hides HITRAN output
            databank="hitran",
            warnings={"AccuracyError": "ignore"},
        )
    except:
        return False

    noise = __loadData(s.get("transmittance_noslit", wunit="nm", Iunit="default"))
    for x in noise:
        noise[x] = noise[x] + np.random.normal(0, 1)

    spectrum = __loadData(s.get("transmittance_noslit", wunit="nm", Iunit="default"))

    # ----- b.) blackbody spectrum of source -----

    spectrum = __sPlanck(spectrum, data["source"])

    # ----- c.) transmission spectrum of windows/beamsplitter -----

    # Beamsplitter
    if data["beamsplitter"] == "AR_ZnSe":
        spectrum = __AR_ZnSe(spectrum)
    elif data["beamsplitter"] == "AR_CaF2":
        spectrum = __AR_CaF2(spectrum)

    # Cell Windows
    if data["cellWindow"] == "CaF2":
        spectrum = __CaF2(spectrum)
        spectrum = __CaF2(spectrum)
    elif data["cellWindow"] == "ZnSe":
        spectrum = __ZnSe(spectrum)
        spectrum = __ZnSe(spectrum)

    # ----- d.) detector response spectrum -----

    if data["detector"] == "MCT":
        spectrum = __ZnSe(spectrum)
        spectrum = __MCT(spectrum)
    elif data["detector"] == "InSb":
        spectrum = __sapphire(spectrum)
        spectrum = __InSb(spectrum)

    return spectrum


if __name__ == "__main__":
    app.run(debug=True)
