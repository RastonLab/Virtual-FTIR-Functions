import radis
from radis import SerialSlabs, Spectrum, calc_spectrum, MergeSlabs
from specutils.fitting import find_lines_threshold

from functions import __sPlanck, __CaF2, __ZnSe, __sapphire, __AR_ZnSe, __AR_CaF2, __InSb, __MCT, __zeroY, __calc_wstep, __multiscan
# ------------------------------
# ----- Spectrum Processing -----
# ------------------------------
def __process_spectrum(params, raw_spectrum):
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
    spec_sPlanck.normalize(normalize_how="mean", inplace=True, force=True)

    # processing for anti-reflective zinc selenide (AR_ZnSe) beamsplitter
    spec_AR_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_ZnSe",
    )
    spec_AR_ZnSe.normalize(normalize_how="mean", inplace=True, force=True)

    # processing for anti-reflective calcium fluoride (AR_CaF2) beamsplitter
    spec_AR_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_CaF2",
    )
    spec_AR_CaF2.normalize(normalize_how="mean", inplace=True, force=True)

    # processing for calcium fluoride (CaF2) cell window
    spec_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="CaF2",
    )
    spec_CaF2.normalize(normalize_how="mean", inplace=True, force=True)


    # processing for zinc selenide (ZnSe) cell window
    spec_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="ZnSe",
    )
    spec_ZnSe.normalize(normalize_how="mean", inplace=True, force=True)


    # processing for sapphire window before detector
    spec_sapphire = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sapphire(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="sapphire",
    )
    spec_sapphire.normalize(normalize_how="mean", inplace=True, force=True)

    # processing for Mercury-Cadmium-Telluride (MCT) detector
    spec_MCT = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __MCT(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="MCT",
    )
    # spec_MCT.normalize(normalize_how="mean", inplace=True, force=True)

    # processing for indium antimonide (InSb) detector
    spec_InSb = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __InSb(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="InSb",
    )
    # spec_InSb.normalize(normalize_how="mean", inplace=True, force=True)

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
            spec_MCT = __multiscan(spec_MCT, params["scan"])
            spec_MCT.normalize(normalize_how="mean", inplace=True, force=True)
            slabs.extend([spec_ZnSe, spec_MCT])
        case "InSb":
            spec_InSb = __multiscan(spec_InSb, params["scan"])
            spec_InSb.normalize(normalize_how="mean", inplace=True, force=True)
            slabs.extend([spec_sapphire, spec_InSb])

    # SerialSlabs() multiplies the transmittance values (y-values) of the selected spectra
    #   https://radis.readthedocs.io/en/latest/source/radis.los.slabs.html#radis.los.slabs.SerialSlabs
    spectrum = SerialSlabs(*slabs, modify_inputs="True")

    # return processed spectrum
    return spectrum


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


def __find_peaks(x_data, y_data, lowerbound, upperbound, threshold=0):
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
