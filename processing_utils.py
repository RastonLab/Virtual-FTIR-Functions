import numpy as np
import warnings

from radis.spectrum.operations import add_array
from radis import Spectrum

from pydantic import ConfigDict, validate_arguments

# filters out specific warning messages
warnings.filterwarnings("ignore", message="invalid value encountered in power")
warnings.filterwarnings("ignore", message="overflow encountered in power")

# -------------------------------------
# ------------- blackbody -------------
# -------------------------------------
# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __sPlanck(spectrum: np.ndarray, source_temp: int) -> np.ndarray:
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
    
    return (0.2 * H * (C**2)) / ((1 / (spectrum * (10**2))) ** 5) * (1 / (np.exp((H * C) / ((1 / (spectrum * (10**2))) * K_B * source_temp)) - 1))


# --------------------------------------
# --------------- window ---------------
# --------------------------------------
# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __CaF2(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a CaF2 cell window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a CaF2 cell window
    """

    return (0.93091) / (1 + (11.12929 / (10000 / spectrum)) ** -12.43933) ** 4.32574


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __ZnSe(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a ZnSe cell window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a ZnSe cell window
    """

    x_um = 10000 / spectrum
    return (0.71015) / ((1 + (20.99353 / x_um) ** -19.31355) ** 1.44348) + -0.13265 / (
        2.25051 * np.sqrt(np.pi / (4 * np.log(2)))
    ) * np.exp(-4 * np.log(2) * ((x_um - 16.75) ** 2) / (2.25051**2))


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __sapphire(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a sapphire window.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a sapphire window
    """

    return 0.78928 / (1 + (11.9544 / (10000 / spectrum)) ** -12.07226) ** 6903.57039


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __AR_ZnSe(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a AR_ZnSe beamsplitter.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a AR_ZnSe beamsplitter
    """

    x_um = 10000 / spectrum
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


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __AR_CaF2(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a AR_CaF2 beamsplitter.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a AR_CaF2 beamsplitter
    """

    x_um = 10000 / spectrum
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
# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __InSb(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for an InSb detector.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with an InSb detector
    """

    x_um = 10000 / spectrum
    return 1.85314E11 * (1 / (1 + np.exp(-(x_um - 5.39001) / 1.80975))) * (
        1 - 1 / (1 + np.exp(-(x_um - 5.39001) / 0.116))
    ) + (3.3E10) / (1.77143 * np.sqrt(np.pi / (4 * np.log(2)))) * np.exp(
        -4 * np.log(2) * ((x_um - 5) ** 2) / (1.77143**2)
    )


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def __MCT(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values for a MCT detector.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a MCT detector
    """

    x_um = 10000 / spectrum # swap spectrum and 10000
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
# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def zeroY(spectrum: np.ndarray) -> np.ndarray:
    """
    Calculates the y-values (y = 1) for background samples.

            Parameters:
                spectrum: An array of x-value for a spectrum

            Returns:
                The y-values associated with a background sample
    """
    return (spectrum * 0) + 1


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def param_check(params: dict[str, object]) -> bool:
    """
    Parses user provided parameters for validity.

        Parameters:
            params (dict): The parameters provided by the user

        Returns:
            True if params are good. Else, returns False
    """

    # check if number of parameters is correct
    if len(params) != 13:
        print("  incorrect amount of params. total params: %s" % (len(params)))
        return False

    # check if parameter names are correct
    valid_params = [
        "beamsplitter",
        "detector",
        "medium",
        "mole",
        "molecule",
        "pressure",
        "resolution",
        "scan",
        "source",
        "waveMax",
        "waveMin",
        "window",
        "zeroFill",
    ]

    for key, value in params.items():
        if (key not in valid_params) or (params[key] is None):
            print(f"  error with key: {key}. Value is: {value}")
            return False

    return True


# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def calc_wstep(resolution: float, zero_fill: int) -> float:
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

        case 0.03125:
            match zero_fill:
                case 0:
                    wstep = 0.01506
                case 1:
                    wstep = 0.00753
                case 2:
                    wstep = 0.003765

        case 0.015625:
            match zero_fill:
                case 0:
                    wstep = 0.00753
                case 1:
                    wstep = 0.003765
                case 2:
                    wstep = 0.001883

    return wstep

# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def multiscan(spectrum: Spectrum, num_scans: int) -> Spectrum:
    '''
    Adds noise to the provided spectrum in chunks to optimize memory usage.

        Parameters:
            spectrum (Spectrum): the spectrum to add noise to
            num_scans (int): the number of scans being run on the sample

        Returns:
            the spectrum with appropriate noise added
    '''
    # add random noise to spectrum
    #   https://radis.readthedocs.io/en/latest/source/radis.spectrum.operations.html#radis.spectrum.operations.add_array
    w = spectrum.get_wavenumber()

    # the maximum scans done per iteration
    scans_per_group = 10
    # how many maximized iterations
    groups = num_scans // scans_per_group

    low = 0
    high = 0.005

    # Adds noise in chunks to save memory space
    for _ in range(groups):
        spectrum = add_array(
            spectrum,
            sum(np.random.normal(low, high, (scans_per_group, len(w)))) / num_scans,
            var="transmittance_noslit",
        )
    
    # does the remaining scans when the number of scans does not evenly divide into 10 (scans_per_group)
    # ex. 115 scans -> the first 110 are done above; the last 5 are done here
    if scans_per_group * groups < num_scans:
        # the number of scans done above == scans_per_group * groups
        # diff == the number of scans remaining
        diff = num_scans - (scans_per_group * groups)
        spectrum = add_array(
            spectrum,
            sum(np.random.normal(low, high, (diff, len(w)))) / num_scans,
            var="transmittance_noslit",
        ) 

    return spectrum

# @validate_arguments(config=ConfigDict(strict=True, arbitrary_types_allowed=True))
def get_component_spectra(w: np.ndarray, source_temp: int) -> tuple[Spectrum, Spectrum, 
                                                                     Spectrum, Spectrum, 
                                                                     Spectrum, Spectrum, 
                                                                     Spectrum, Spectrum]:
    '''
    Calculates the spectra for the components of the spectrometer.

        Parameters:
            w (np.ndarray): the x-values for all of the spectra
            source_temp (int): the source temperature for the blackbody spectrum

        Returns:
            a tuple containing all of the component spectra
    '''
    # processing for blackbody spectrum (sPlanck)
    spec_sPlanck = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sPlanck(w, source_temp)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="sPlanck",
    )
    # Normalize the blackbody spectrum to 1
    spec_sPlanck.normalize(normalize_how="max", inplace=True, force=True)

    # processing for anti-reflective zinc selenide (AR_ZnSe) beamsplitter
    spec_AR_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_ZnSe(w) ** (3/5)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="AR_ZnSe",
    )

    # processing for anti-reflective calcium fluoride (AR_CaF2) beamsplitter
    spec_AR_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_CaF2(w) ** (3/5)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="AR_CaF2",
    )

    # processing for calcium fluoride (CaF2) cell window
    spec_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __CaF2(w) ** (2/5)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="CaF2",
    )

    # processing for zinc selenide (ZnSe) cell window
    spec_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __ZnSe(w) ** (2/5)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="ZnSe",
    )

    # processing for sapphire window before detector
    spec_sapphire = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sapphire(w) ** (1/5)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="sapphire",
    )

    # processing for Mercury-Cadmium-Telluride (MCT) detector
    spec_MCT = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __MCT(w)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="MCT",
    )
    # Normalize the MCT spectrum to 1
    spec_MCT.normalize(normalize_how="max", inplace=True, force=True)

    # processing for indium antimonide (InSb) detector
    spec_InSb = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __InSb(w)},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="InSb",
    )
    # Normalize the InSb spectrum to 2
    spec_InSb.normalize(normalize_how="max", inplace=True, force=True)

    _, y_value = spec_InSb.get("transmittance_noslit")
    y_value *= 2

    spec_InSb = Spectrum(
        {"wavenumber": w, "transmittance_noslit": y_value},
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="InSb",
    )
    # End InSb Normalization
    
    return (spec_sPlanck, spec_AR_ZnSe, spec_AR_CaF2, spec_CaF2, spec_ZnSe, 
            spec_sapphire, spec_MCT, spec_InSb)