import json
import numpy as np
import matplotlib.pyplot as plt

from functions import __param_check, __generate_spectra


def main():
    # read local data file into a dictionary
    with open("../data.json", "r") as data_json:
        data = json.load(data_json)

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

    xs = []
    ys = []
    for key in result:
        xs.append(float(key))
        ys.append(float(result[key]))
    plt.plot(np.array(xs), np.array(ys), "blue")
    plt.show()
    print("----- finished plot, ending -----")


if __name__ == "__main__":
    main()
