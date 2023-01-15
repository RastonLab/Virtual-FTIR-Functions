import astropy.units as u
import numpy as np
import radis
from radis import SerialSlabs, Spectrum, calc_spectrum
from radis.spectrum.operations import add_array, multiply
from specutils import SpectralRegion
from specutils.fitting import find_lines_derivative, find_lines_threshold
from specutils.manipulation import noise_region_uncertainty


# -------------------------------------
# ------------- blackbody -------------
# -------------------------------------
def __sPlanck(spectrum, source_temp):
    """
    Calculates the y-values of a Blackbody spectrum.

            Parameters:
                spectrum: An array of x-value for a spectrum
                source_temp: The source temperature associated with the spectrum

            Returns:
                The y-values of a Blackbody spectrum
    """

    H = 6.62606957e-34
    C = 2.99792458e8
    K_B = 1.3806488e-23

    return ((0.2 * H * (C**2)) / (((spectrum * (10**-9)) ** 4) * spectrum)) * (
        1 / (np.exp((H * C) / ((spectrum * (10**-9)) * K_B * source_temp)) - 1)
    )


# --------------------------------------
# --------------- window ---------------
# --------------------------------------
def __CaF2(spectrum):
    """
    Calculates the y-values for a CaF2 cell window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a CaF2 cell window
    """

    return (0.93091) / (1 + (11.12929 / (spectrum / 1000)) ** -12.43933) ** 4.32574


def __ZnSe(spectrum):
    """
    Calculates the y-values for a ZnSe cell window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a ZnSe cell window
    """

    x_um = spectrum / 1000
    return (0.71015) / ((1 + (20.99353 / x_um) ** -19.31355) ** 1.44348) + -0.13265 / (
        2.25051 * np.sqrt(np.pi / (4 * np.log(2)))
    ) * np.exp(-4 * np.log(2) * ((x_um - 16.75) ** 2) / (2.25051**2))


def __sapphire(spectrum):
    """
    Calculates the y-values for a sapphire window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a sapphire window
    """

    return 0.78928 / (1 + (11.9544 / (spectrum / 1000)) ** -12.07226) ** 6903.57039


def __AR_ZnSe(spectrum):
    """
    Calculates the y-values for a AR_ZnSe beamsplitter.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a AR_ZnSe beamsplitter
    """

    x_um = spectrum / 1000
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


def __AR_CaF2(spectrum):
    """
    Calculates the y-values for a AR_CaF2 beamsplitter.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a AR_CaF2 beamsplitter
    """

    x_um = spectrum / 1000
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


# --------------------------------------
# -------------- detector --------------
# --------------------------------------
def __InSb(spectrum):
    """
    Calculates the y-values for an InSb detector.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with an InSb detector
    """

    x_um = spectrum / 1000
    return 1.97163e11 * (1 / (1 + np.exp(-(x_um - 5.3939) / 1.6624))) * (
        1 - 1 / (1 + np.exp(-(x_um - 5.3939) / 0.11925))
    ) + (3.3e10) / (2.44977 * np.sqrt(np.pi / (4 * np.log(2)))) * np.exp(
        -4 * np.log(2) * ((x_um - 5) ** 2) / (2.44977**2)
    )


def __MCT(spectrum):
    """
    Calculates the y-values for a MCT detector.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a MCT detector
    """

    x_um = spectrum / 1000
    return (
        (1.98748 * (10**9))
        + (2.10252 * (10**10))
        * (1 / (1 + np.exp(-(x_um - 20.15819) / 5.73688)))
        * (1 - 1 / (1 + np.exp(-(x_um - 20.15819) / 1.11659)))
        + (1.3 * (10**9))
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 18.6) ** 2) / (2**2))
    )


# -------------------------------------
# ---------- helper functions ----------
# ------------------------------------
def __zeroY(spectrum):
    """
    Calculates the y-values (y = 1) for background samples.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a background sample
    """
    return (spectrum * 0) + 1


def __param_check(params):
    """
    Parses user provided parameters for validity.

        Parameters:
            params (dict): The parameters provided by the user

        Returns:
            True if params are good. Else, returns False
    """

    # check if number of parameters is correct
    if len(params) != 11:
        print("  not enough params. total params: %s" % (len(params)))
        return False

    # check if parameter names are correct
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

    for key, value in params.items():
        if (key not in valid_params) or (params[key] is None):
            print(f"  error with key: {key}. Value is: {value}")
            return False

    return True


def __calc_wstep(resolution, zero_fill):
    """
    Calculates the appropriate wstep for a spectrum based on the given resolution and zero fill.

        Parameters:
            resolution (int): the given resolution
            zero_fill (int): the given zero fill

        Returns:
            The calculated wstep
    """

    match resolution:
        case 1:
            match zero_fill:
                case 0:
                    wstep = 0.481927711
                case 1:
                    wstep = 0.240963855
                case 2:
                    wstep = 0.120481928

        case 0.5:
            match zero_fill:
                case 0:
                    wstep = 0.240963855

                case 1:
                    wstep = 0.120481928

                case 2:
                    wstep = 0.060240964

        case 0.25:
            match zero_fill:
                case 0:
                    wstep = 0.120481928
                case 1:
                    wstep = 0.060240964
                case 2:
                    wstep = 0.030120482

        case 0.125:
            match zero_fill:
                case 0:
                    wstep = 0.060240964
                case 1:
                    wstep = 0.030120482
                case 2:
                    wstep = 0.015060241

        case 0.0625:
            match zero_fill:
                case 0:
                    wstep = 0.030120482
                case 1:
                    wstep = 0.015060241
                case 2:
                    wstep = 0.00753012

    return wstep


# ------------------------------
# ----- Spectrum Processing -----
# ------------------------------
def __process_spectrum(params, raw_spectrum, find_peaks):
    """
    The following function takes a 'raw spectrum' generated using Radis's
    'calc_spectrum()' function and performing custom equations that virtualize
    the spectrum in a spectrometer.

    This is achieved by creating additional spectra based on functions that
    approximate the behavior of FTIR components. Those spectra are then
    multiplied into the base spectrum.

    All steps except pre-processing, Step A, and post-processing are repeated
    a number of times given by the user. At the end of each loop, the spectrum
    is normalized.

        Parameters:
            params (dict): The parameters provided by the user
            raw_spectrum (Spectrum object): The spectrum generated from 'calc_spectrum()'
            find_peaks (boolean): Tells wether or not the find peaks algorithm should run

        Returns:
            The processed spectrum as a dictionary

    """

    # ----- Pre-processing -----
    # generate the necessary spectra for blackbody, beamsplitters, cell windows, detectors. the spectra are generated based on the function provided in the call to the Spectrum constructor

    # returns the x-values of calc_spectrum() in an array
    #   https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.get_wavenumber
    w = raw_spectrum.get_wavenumber()

    # processing for blackbody spectrum (sPlanck)
    spec_sPlanck = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sPlanck(w, params["source"])},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="sPlanck",
    )

    # processing for anti-reflective zinc selenide (AR_ZnSe) beamsplitter
    spec_AR_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_ZnSe",
    )

    # processing for anti-reflective calcium fluoride (AR_CaF2) beamsplitter
    spec_AR_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_CaF2",
    )

    # processing for calcium fluoride (CaF2) cell window
    spec_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="CaF2",
    )

    # processing for zinc selenide (ZnSe) cell window
    spec_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="ZnSe",
    )

    # processing for sapphire window before detector
    spec_sapphire = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sapphire(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="sapphire",
    )

    # processing for Mercury-Cadmium-Telluride (MCT) detector
    spec_MCT = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __MCT(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="MCT",
    )

    # processing for indium antimonide (InSb) detector
    spec_InSb = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __InSb(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="InSb",
    )

    # list of spectra to multiply
    slabs = []

    # ----- a.) transmission spectrum of gas sample -----
    slabs.append(raw_spectrum)

    # ----- b.) blackbody spectrum of source -----
    slabs.append(spec_sPlanck)

    # ----- c.) transmission spectrum of windows/beamsplitter -----
    # ----- c.1) Beamsplitter -----
    match params["beamsplitter"]:
        case "AR_ZnSe":
            slabs.append(spec_AR_ZnSe)
        case "AR_CaF2":
            slabs.append(spec_AR_CaF2)

    # ----- c.2) cell windows -----
    match params["cellWindow"]:
        case "CaF2":
            slabs.extend([spec_CaF2, spec_CaF2])
        case "ZnSe":
            slabs.extend([spec_ZnSe, spec_ZnSe])

    # ----- d.) detector response spectrum -----
    match params["detector"]:
        case "MCT":
            slabs.extend([spec_ZnSe, spec_MCT])
        case "InSb":
            slabs.extend([spec_sapphire, spec_InSb])

    # SerialSlabs() multiplies the transmittance values (y-values) of the selected spectra
    #   https://radis.readthedocs.io/en/latest/source/radis.los.slabs.html#radis.los.slabs.SerialSlabs
    spectrum = SerialSlabs(*slabs, modify_inputs="True")

    spectrum.normalize(normalize_how="mean", inplace=True, force=True)

    # add random noise to spectrum
    #   https://radis.readthedocs.io/en/latest/source/radis.spectrum.operations.html#radis.spectrum.operations.add_array
    spectrum = add_array(
        spectrum,
        sum(np.random.normal(0, 1, (params["numScan"], len(w))))
        / params["numScan"],
        var="transmittance_noslit",
    )

    # Post-processing - Find Peaks
    # Not done on background samples
    # https://radis.readthedocs.io/en/latest/auto_examples/plot_specutils_processing.html#sphx-glr-auto-examples-plot-specutils-processing-py
    # if find_peaks:
    #     find_peaks = spectrum.to_specutils()
    #     noise_region = SpectralRegion(
    #         (1 / data["minWave"]) / u.cm, (1 / data["maxWave"]) / u.cm
    #     )
    #     find_peaks = noise_region_uncertainty(find_peaks, noise_region)
    #     lines = find_lines_threshold(find_peaks, noise_factor=6)
    #     print()
    #     print(lines)

    # return processed spectrum
    return spectrum


def __generate_spectrum(params):
    """
    Generates a spectrum using Radis's 'calc_spectrum()' function based
    on user parameters. That spectrum is then processed by
    '__process_spectrum()'.

    If there is an issue with the Radis library, the error message is returned.

        Parameters:
            params (dict): The parameters provided by the user

        Return:
            The raw spectrum, or the message text if an error occurs
    """

    # resolution of wavenumber grid (cm^-1)
    #   https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum:~:text=wstep%20(float%20(,%27auto%27)
    wstep = __calc_wstep(params["resolution"], params["zeroFill"])

    try:
        # ----- a.) transmission spectrum of gas sample -----
        #   https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum
        spectrum = calc_spectrum(
            params["minWave"],
            params["maxWave"],
            molecule=params["molecule"],
            isotope="1,2,3",
            pressure=params["pressure"],
            Tgas=294.15,
            path_length=10,
            wstep=wstep,
            databank="hitran",
            verbose=False,
            warnings={"AccuracyError": "ignore"},
        )
    except radis.misc.warning.EmptyDatabaseError:
        return None, True, "error: No line in the specified wavenumber range"
    except Exception as e:
        match str(e):
            case "Failed to retrieve data for given parameters.":
                return (
                    None,
                    True,
                    "error: HITRAN data does not exist for requested molecule.",
                )
            case other:
                return None, True, str(e)

    return spectrum, False, None


def __process_background(raw_spectrum):
    """
    Accepts a spectrum generated using '__generate_spectrum()'.
    A background by default has all y-values of one.

        Parameters:
            raw_spectrum (Spectrum object): The spectrum generated from 'calc_spectrum()'

        Return:
            The processed background sample with y-values of one
    """

    spec_zeroY = Spectrum(
        {
            "wavenumber": raw_spectrum.get_wavenumber(),
            "transmittance_noslit": __zeroY(raw_spectrum.get_wavenumber()),
        },
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="Background",
    )

    return spec_zeroY
