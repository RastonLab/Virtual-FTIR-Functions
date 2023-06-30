import radis
from radis import SerialSlabs, Spectrum, calc_spectrum, MergeSlabs
from specutils.fitting import find_lines_threshold
from functions import zeroY, calc_wstep, multiscan, get_component_spectra
# ------------------------------
# ----- Spectrum Processing -----
# ------------------------------
def process_spectrum(params, raw_spectrum):
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

    spec_sPlanck, spec_AR_ZnSe, spec_AR_CaF2, spec_CaF2, spec_ZnSe, spec_sapphire, \
        spec_MCT, spec_InSb = get_component_spectra(w, params["source"])

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
    match params["window"]:
        case "CaF2":
            slabs.extend([spec_CaF2, spec_CaF2])
        case "ZnSe":
            slabs.extend([spec_ZnSe, spec_ZnSe])

    # ----- d.) detector response spectrum -----
    match params["detector"]:
        case "MCT":
            # spec_MCT = multiscan(spec_MCT, params["scan"])
            slabs.extend([spec_ZnSe, spec_MCT])
        case "InSb":
            # spec_InSb = multiscan(spec_InSb, params["scan"])
            slabs.extend([spec_sapphire, spec_InSb])

    # SerialSlabs() multiplies the transmittance values (y-values) of the selected spectra
    #   https://radis.readthedocs.io/en/latest/source/radis.los.slabs.html#radis.los.slabs.SerialSlabs
    spectrum = SerialSlabs(*slabs, modify_inputs="True")
    # spectrum = multiscan(spectrum, params["scan"])

    # return processed spectrum
    return spectrum


def process_background(raw_spectrum):
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
            "transmittance_noslit": zeroY(raw_spectrum.get_wavenumber()),
        },
        wunit="cm-1",
        units={"transmittance_noslit": ""},
        name="Background",
    )

    return spec_zeroY


def generate_spectrum(params):
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
    wstep = calc_wstep(params["resolution"], params["zeroFill"])

    try:
        # ----- a.) transmission spectrum of gas sample -----
        #   https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum
        spectrum = calc_spectrum(
            params["waveMin"],
            params["waveMax"],
            molecule=params["molecule"],
            isotope="1,2,3",
            pressure=params["pressure"],
            Tgas=294.15,
            path_length=10,
            wstep=wstep,
            databank="hitran",
            verbose=False,
            warnings={
                "AccuracyError": "ignore",
                "AccuracyWarning": "ignore"},
            mole_fraction={params["molecule"]: params["mole"]},
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


def find_peaks(x_data, y_data, lowerbound, upperbound, threshold=0):
    try:
        spectrum = Spectrum.from_array(
            x_data, y_data, "absorbance_noslit", wunit="cm-1", unit=""
        )
        new_spec = (
            spectrum.to_specutils()
        )
        lines = find_lines_threshold(new_spec, noise_factor=1)
    except:
        return None

    peaks = {}
    for num, peak_type, _ in lines:
        index = x_data.index(float(num.value))
        if x_data[index] >= lowerbound and x_data[index] <= upperbound:
            if peak_type == "emission" and y_data[index] >= threshold:
                peaks[round(float(num.value), 4)] = round(y_data[index], 4)

    return peaks
