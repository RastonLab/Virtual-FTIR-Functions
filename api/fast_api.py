# FastAPI server
#   https://fastapi.tiangolo.com/
#   uvicorn fast_api:app --reload

# https://realpython.com/fastapi-python-web-apis/
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Literal

from functions import __param_check, __generate_spectra, __generate_background

molecule_list = Literal[
    "C2H2",
    "C2H4",
    "C2H6",
    "C2N2",
    "C4H2",
    "CF4",
    "CH3Br",
    "CH3Cl",
    "CH3CN",
    "CH3OH",
    "CH4",
    "ClO",
    "ClONO2",
    "CO",
    "CO2",
    "COCl2",
    "COF2",
    "CS",
    "H2",
    "H2CO",
    "H2O",
    "H2O2",
    "H2S",
    "HBr",
    "HC3N",
    "HCl",
    "HCN",
    "HCOOH",
    "HF",
    "HI",
    "HNO3",
    "HO2",
    "HOBr",
    "HOCl",
    "N2",
    "N2O",
    "NH3",
    "NO",
    "NO+",
    "NO2",
    "O",
    "O2",
    "O3",
    "OCS",
    "OH",
    "PH3",
    "SF6",
    "SO2",
    "SO3",
]


class Payload(BaseModel):
    minWave: float
    maxWave: float
    molecule: molecule_list
    pressure: float
    resolution: float
    numScan: int
    zeroFill: Literal[0, 1, 2]
    source: Literal[1700, 3100]
    beamsplitter: Literal["AR_ZnSe", "AR_CaF2"]
    cellWindow: Literal["ZnSe", "CaF2"]
    detector: Literal["MCT", "InSb"]


app = FastAPI()
origins = [
    "http://localhost:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    return "<h1>ftir</h1>"


@app.post("/fetch_background/")
async def create_item(payload: Payload):
    data = {
        "minWave": payload.minWave,
        "maxWave": payload.maxWave,
        "molecule": payload.molecule,
        "pressure": payload.pressure,
        "resolution": payload.resolution,
        "numScan": payload.numScan,
        "zeroFill": payload.zeroFill,
        "source": payload.source,
        "beamsplitter": payload.beamsplitter,
        "cellWindow": payload.cellWindow,
        "detector": payload.detector,
    }

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


@app.post("/post_json/")
async def create_item(payload: Payload):
    data = {
        "minWave": payload.minWave,
        "maxWave": payload.maxWave,
        "molecule": payload.molecule,
        "pressure": payload.pressure,
        "resolution": payload.resolution,
        "numScan": payload.numScan,
        "zeroFill": payload.zeroFill,
        "source": payload.source,
        "beamsplitter": payload.beamsplitter,
        "cellWindow": payload.cellWindow,
        "detector": payload.detector,
    }

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
