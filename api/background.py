import json
import numpy as np
import matplotlib.pyplot as plt

from functions import __param_check, __generate_background


def main():
    # read local data file into a dictionary
    with open("./data.json", "r") as data_json:
        data = json.load(data_json)

    print(data)

    # verify the information in the dictionary
    __param_check(data)

    print("----- output verified params to console as self-check -----")
    for key, value in data.items():
        print("  %s: %s" % (key, value))

    # generate background noise that will be added to calc_spectrum result
    print("----- start __generate_background() -----")
    result = __generate_background(data)
    print("----- end __generate_background() -----")

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
