from radis import calc_spectrum, Spectrum, SerialSlabs
from radis.spectrum.operations import add_array, multiply

from specutils import SpectralRegion
from specutils.manipulation import noise_region_uncertainty
from specutils.fitting import find_lines_threshold, find_lines_derivative

import astropy.units as u
import numpy as np

# --------------------------------------
# ---- spectra calculation functions ----
# --------------------------------------

def __CaF2(data):
    '''
    Returns a function that approximates a CaF2 Cell Window.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    return (0.93091) / (1 + (11.12929 / (data / 1000)) ** -12.43933) ** 4.32574


def __ZnSe(data):
    '''
    Returns a function that approximates a ZnSe Cell Window.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    x_um = data / 1000
    return (0.71015) / ((1 + (20.99353 / x_um) ** -19.31355) ** 1.44348) + -0.13265 / (
        2.25051 * np.sqrt(np.pi / (4 * np.log(2)))
    ) * np.exp(-4 * np.log(2) * ((x_um - 16.75) ** 2) / (2.25051**2))


def __sapphire(data):
    '''
    Returns a function that approximates a Sapphire Window.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    return 0.78928 / (1 + (11.9544 / (data / 1000)) ** -12.07226) ** 6903.57039


def __AR_ZnSe(data):
    '''
    Returns a function that approximates a AR_ZnSe Beamsplitter.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    x_um = data / 1000
    return (
        (0.82609) / ((1 + ((34.63971 / x_um) ** -8.56269)) ** 186.34792)
        + -0.47
        / (0.55 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 1.47) ** 2) / (0.55**2))
        + -0.03456
        / (0.4 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.88) ** 2) / (0.4**2))
        + -0.009
        / (0.3 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 6.16) ** 2) / (0.3**2))
        + -0.09
        / (1 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 16.2) ** 2) / (1**2))
        + -0.08
        / (1 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 17.4) ** 2) / (1**2))
        + 1.12
        / (8 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 9.5) ** 2) / (8**2))
        + 0.11546
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 4.9) ** 2) / (2**2))
        + 0.21751
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.6) ** 2) / (2**2))
        + -0.05
        / (0.07 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 0.8) ** 2) / (0.07**2))
    )


def __AR_CaF2(data):
    '''
    Returns a function that approximates a AR_CaF2 Beamsplitter.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    x_um = data / 1000
    return (
        (0.9795) / ((1 + ((18.77617 / x_um) ** -6.94246)) ** 91.98745)
        + -0.06
        / (0.08 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 0.76) ** 2) / (0.08**2))
        + -0.06
        / (0.2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * (x_um - 1.06) ** 2 / 0.20**2)
        + -0.6
        / (3.0 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 4.85) ** 2) / (3.0**2))
        + -0.35
        / (1.0 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 9.40) ** 2) / (1.00**2))
        + 0.05
        / (0.8 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.60) ** 2) / (0.8**2))
        + 0.04
        / (0.5 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 7.75) ** 2) / (0.50**2))
        + -0.01
        / (0.6 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 6.55) ** 2) / (0.6**2))
        + 0.01
        / (0.5 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 1.82) ** 2) / (0.5**2))
    )


def __InSb(data):
    '''
    Returns a function that approximates a InSb Dectector.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    x_um = data / 1000
    return 1.97163e11 * (1 / (1 + np.exp(-(x_um - 5.3939) / 1.6624))) * (
        1 - 1 / (1 + np.exp(-(x_um - 5.3939) / 0.11925))
    ) + (3.3e10) / (2.44977 * np.sqrt(np.pi / (4 * np.log(2)))) * np.exp(
        -4 * np.log(2) * ((x_um - 5) ** 2) / (2.44977**2)
    )


def __MCT(data):
    '''
    Returns a function that approximates a MCT Dectector.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''

    x_um = data / 1000
    return (
        (1.98748 * (10**9))
        + (2.10252 * (10**10))
        * (1 / (1 + np.exp(-(x_um - 20.15819) / 5.73688)))
        * (1 - 1 / (1 + np.exp(-(x_um - 20.15819) / 1.11659)))
        + (1.3 * (10**9))
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 18.6) ** 2) / (2**2))
    )

def __zeroY(data):
    '''
    Returns a function of y = 1 for generating background samples.

            Parameters:
                data (int): the range of x-values for the function
            
            Returns:
                The approximated function
    '''
    return (data * 0) + 1

def __sPlanck(spectrum, temp):
    '''
    Returns a function that generates a Blackbody Spectrum.

            Parameters:
                data (int): the range of x-values for the spectrum
            
            Returns:
                The Blackbosy Spectrum
    '''

    H = 6.62606957e-34
    C = 2.99792458e8
    K_B = 1.3806488e-23

    return ((0.2 * H * (C**2)) / (((spectrum * (10**-9)) ** 4) * spectrum)) * (
        1 / (np.exp((H * C) / ((spectrum * (10**-9)) * K_B * temp)) - 1)
    )


# -------------------------------------
# ---------- helper functions ----------
# ------------------------------------
def __loadData(s):
    '''
    Takes a spectrum and converts its data to a dictionary of x and y values.

            Parameters:
                s (Spectrum): The spectrum to convert
            
            Returns:
                A dictionary of x (key) and y (value) values
    '''

    data = {}

    for key, val in zip(s[0], s[1]):
        data[float(key)] = float(val)

    return data


def __param_check(data):
    '''
    Checks user given parameters and makes sure they are vaild.

        Parameters:
            data (dict): the user given parameters

        Returns:
            Good if params are good. Else, returns an error message
    '''

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
    '''
    Calculates the apropriate wstep for a spectrum based on the given resolution and zero fill.

        Parameters:
            resolution (int): the given resolution
            xero_fill (int): the given zero fill

        Returns:
            The calculated wstep
    '''

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

# ------------------------------
# ----- Spectra Processing -----
# ------------------------------

def __find_peaks(data, spectrum):    
    find_peaks = spectrum.to_specutils('transmittance_noslit', 'nm_vac') #had to put vac, cannot put air
    lines = find_lines_threshold(find_peaks, noise_factor=200000000) #what should this be set to?

    lines.write ("lines.txt", format="ascii", overwrite=True)


def __process_spectra(data, s, find_peaks):
    '''
    Takes a spectrum and adds noise that approximates a real spectrometer.

    This is achieved by creating additonal spectra based on functions that 
    approximate the behavior of FTIR components. Those spectra are then 
    multiplied into the base spectrum. 

    All steps except pre-processing, Step A, and post-processing are repeated
    a number of times given by the user. At the end of each loop, the spectrum
    is normalized.

    If the spectrum is a background sample, no post-processing is done. 
    Otherwise, an algorithm for find the peaks of the spectrum is executed

        Parameters:
            data (dict): User given parameters
            s (Spectrum): The given spectrum
            find_peaks (boolean): Tells wether or not the find peaks algorithm should run

        Returns:
            The processed spectrum as a dictionary

    '''

    # ----- Pre-processing -----
    # Generate the necessary spectra for each component of the following processing steps
    # The spectra are generated based on the function provided in the call to the Spectrum constructor
    w = s.get_wavelength()

    # TODO The Radis dev gave us this part
    spec_sPlanck = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __sPlanck(w, data["source"]),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_AR_ZnSe = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __AR_ZnSe(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_AR_CaF2 = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __AR_CaF2(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_CaF2 = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __CaF2(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_ZnSe = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __ZnSe(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_sapphire = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __sapphire(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_MCT = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __MCT(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    spec_InSb = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __InSb(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )
    
    # ----- b.) blackbody spectrum of source -----
    # SerialSlabs() multiplies the transmittance values (y-values) of the provided spectrum, s, with 
    # the corresponding transmittance values of the spec_sPlanck spectrum
    # https://radis.readthedocs.io/en/latest/source/radis.los.slabs.html#radis.los.slabs.SerialSlabs
    spectrum = SerialSlabs(s, spec_sPlanck)

    # This loop simpulates scans and runs as many times as the user indicated in "Number of Scans"
    for x in range(data["numScan"]):
        
        # ----- c.) transmission spectrum of windows/beamsplitter -----

        # ----- c.1) Beamsplitter -----
        if data["beamsplitter"] == "AR_ZnSe":
            spectrum = SerialSlabs(spectrum, spec_AR_ZnSe)
        elif data["beamsplitter"] == "AR_CaF2":
            spectrum = SerialSlabs(spectrum, spec_AR_CaF2)

        # ----- c.2) Cell Windows -----
        if data["cellWindow"] == "CaF2":
            spectrum = SerialSlabs(spectrum, spec_CaF2)
            spectrum = SerialSlabs(spectrum, spec_CaF2)
        elif data["cellWindow"] == "ZnSe":
            spectrum = SerialSlabs(spectrum, spec_ZnSe)
            spectrum = SerialSlabs(spectrum, spec_ZnSe)

         # ----- d.) detector response spectrum -----
        if data["detector"] == "MCT":
            spectrum = SerialSlabs(spectrum, spec_ZnSe)
            spectrum = SerialSlabs(spectrum, spec_MCT)
        elif data["detector"] == "InSb":
            spectrum = SerialSlabs(spectrum, spec_sapphire)
            spectrum = SerialSlabs(spectrum, spec_InSb)

        spectrum = add_array(
            spectrum,
            np.random.normal(0, 200000000, len(spectrum)),
            var="transmittance_noslit",
        )

        # Normalize data
        # numbers = __loadData(
        #     spectrum.get("transmittance_noslit", wunit="nm", Iunit="default")
        # )
        # factor = 1 / sum(numbers.values())

        # spectrum = multiply(spectrum, factor, var="transmittance_noslit")

    # Post-processing - Find Peaks
    # Not done on background samples
    # https://radis.readthedocs.io/en/latest/auto_examples/plot_specutils_processing.html#sphx-glr-auto-examples-plot-specutils-processing-py
    if find_peaks:
        __find_peaks(data, spectrum)

    # Return spectrum as a dictionary
    return __loadData(spectrum.get("transmittance_noslit", wunit="nm", Iunit="default"))

def __generate_spectra(data):
    '''
    Generates a spectrum using radis based on user given parameters.

    That spectrum is then processed by __process_spectra.

    If there is an issue reaching out to radis, the error is caught and False is returned

        Paramters:
            data (dict): the user given parameters

        Return:
            The fully processed spectrum or False if there is an error
    '''

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

    return __process_spectra(data, s, find_peaks=True)

def __generate_background(data):
    '''
    Generates a background sample using radis based on user given parameters.

    That sample is then processed by __process_spectra.

    If there is an issue reaching out to radis, the error is caught and False is returned

        Paramters:
            data (dict): the user given parameters

        Return:
            The fully processed sample or False if there is an error
    '''

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

    w = s.get_wavelength()

    spec_zeroY = Spectrum(
        {
            "wavelength": w,
            "transmittance_noslit": __zeroY(w),
            "radiance_noslit": np.zeros_like(w),
        },
        wunit="nm",
        units={"radiance_noslit": "mW/cm2/sr/nm", "transmittance_noslit": ""},
        name="CaF2 window",
    )

    return __process_spectra(data, spec_zeroY, find_peaks=False)