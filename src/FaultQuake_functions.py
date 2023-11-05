# FaultQuake_functions.py

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import os


# DATA = defaultdict(dict)

import os

def export_faults_to_xml(faults, outputfilename):
    # Check if './output_files/Fault_OQ' folder exists, if not, create it
    os.makedirs(f'./{outputfilename}/Fault_OQ', exist_ok=True)

    for fault_name in faults:
        # ScR = faults[fault_name]['ScR']
        # Length = np.append(Length, faul
        xml_filename = f'./{outputfilename}/Fault_OQ/{fault_name}.xml'
        current_directory = os.getcwd()
        print(f'Current directory: {current_directory}')

        # Determine XML filename from fault name and save in './Sources' folder
        xml_filename = f'./{outputfilename}/Fault_OQ/{fault_name}.xml'
        with open(xml_filename, 'w') as fid:
            # Write the header information
            fid.write('<?xml version="1.0" encoding="utf-8"?>\n')
            fid.write('<nrml xmlns="http://openquake.org/xmlns/nrml/0.4" xmlns:gml="http://www.opengis.net/gml">\n')
            fid.write(f'    <sourceModel name="{fault_name}">\n')

            # Write the fault details to the XML
            fid.write(
                f'        <simpleFaultSource id="{faults[fault_name]["id"]}" name="Simple Fault Source" tectonicRegion="Active Shallow Crust">\n')
            fid.write('            <simpleFaultGeometry>\n')
            fid.write('                <gml:LineString>\n')
            fid.write('                    <gml:posList>\n')
            fid.write('                        <!-- Fill in coordinates here -->\n')
            fid.write('                    </gml:posList>\n')
            fid.write('                </gml:LineString>\n')
            fid.write(f'                <dip>{faults[fault_name]["Dip"]}</dip>\n')
            fid.write(
                f'                <upperSeismoDepth>{faults[fault_name]["Seismogenic_Thickness"] - faults[fault_name]["Telap"]:e}</upperSeismoDepth>\n')
            fid.write(f'                <lowerSeismoDepth>{faults[fault_name]["Seismogenic_Thickness"]:e}</lowerSeismoDepth>\n')
            fid.write('            </simpleFaultGeometry>\n')
            fid.write(f'            <magScaleRel>{faults[fault_name]["ScR"]}</magScaleRel>\n')
            fid.write('            <ruptAspectRatio>2.0000000E+00</ruptAspectRatio>\n')  # Placeholder
            # Uncomment the next line if incrementalMFD needs to be included
            # fid.write(f'            <incrementalMFD minMag="{fault["Mmin"]}" binWidth="{fault["bin"]}">\n')
            fid.write(f'                <occurRates>{" ".join(f"{rate:e}" for rate in faults[fault_name]["rates"])}</occurRates>\n')
            fid.write('            </incrementalMFD>\n')
            fid.write('            <rake>9.0000000E+01</rake>\n')  # Placeholder
            fid.write('        </simpleFaultSource>\n')

            # Finish the XML
            fid.write('    </sourceModel>\n')
            fid.write('</nrml>\n')


import numpy as np
import matplotlib.pyplot as plt
import os


def TruncatedGR(fauls, c, d, outputfilename, faultnames, mags, mts, Morates, ids, nfault, bin, bs, DATA):
    outputname = f"{outputfilename}_AR_TruncatedGR.txt"
    with open(os.path.join('./output_files', outputname), 'w') as fidout:
        # print a title, followed by a blank line
        fidout.write('id Mmin bin rates name\n')
        for i in range(nfault):  # cycle for number of faults
            magnitude_range = np.arange(mts[i], mags[i] + bin, bin)
            M = 10 ** (c * magnitude_range + d)
            Beta = (2 / 3) * bs[i]
            Mt = 10 ** (c * mts[i] + d)
            Mxp = 10 ** (c * (mags[i] + bin) + d)  # Assume max magnitude is reached with an additional bin.
            TruncGR = (((Mt / M) ** Beta - (Mt / Mxp) ** Beta) / (1 - (Mt / Mxp) ** Beta))
            Incremental = np.concatenate((np.diff(TruncGR[::-1])[::-1], [TruncGR[-1]]))
            Incremental_Morate = Incremental * M
            Incremental_Morate_balanced = Incremental_Morate * Morates[i] / np.sum(Incremental_Morate)
            cons_tassi_ind = Incremental * Incremental_Morate_balanced / Incremental_Morate
            cumulative_rates = np.cumsum(cons_tassi_ind[::-1])[::-1]
            Mbalanced = np.sum(cons_tassi_ind * M)
            out_Rates = [ids[i], mts[i], bin] + cons_tassi_ind.tolist()

            fault_data = {}
            fault_data = {
                "rates": out_Rates # Your 'i' index
            }

            # Adding the outputs of the moment budget to the faults dictionary
            faults[faultname[i]].update(fault_data)

            # Update DATA dictionary
            DATA[faultnames[i]] = {'id': ids[i], 'rates': cons_tassi_ind, 'bin': bin}

            # Write to output file
            fidout.write(f"{ids[i]}, {mts[i]:3.1f}, {bin:3.1f}, ")
            fidout.write(' '.join(f"{rate:5.4e}" for rate in cons_tassi_ind))
            fidout.write(f", {faultnames[i]}\n")

            # Plotting
            plt.figure(i)
            plt.semilogy(magnitude_range, cumulative_rates, 'ok')
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultnames[i])

            # Saving the figure
            figname = f"{outputfilename}_AR_TruncatedGR_rates_{faultnames[i]}.eps"
            plt.savefig(os.path.join('./output_files', figname), format='epsc')
    export_faults_to_xml(faults, outputfilename)

    aaa=111

    return outputname




def CHGaussPoiss(faults, c, d, outputfilename, faultname, mag, sdmag, Morate, id, nfault, w, Hpois, bin):
    outputname = f"{outputfilename}_AR_ChGaussPoisson_rates.txt"
    outputnameProbability = f"{outputfilename}_AR_ChGaussPoisson_Probability.txt"

    # Create output_files directory if it doesn't exist
    os.makedirs('./output_files/', exist_ok=True)

    # Open files for writing the outputs
    with open(f'./output_files/{outputname}', 'w') as fidout, open(f'./output_files/{outputnameProbability}',
                                                                   'w') as fidoutProb:
        # Print a title, followed by a blank line
        fidout.write('id Mmin bin rates name\n')
        fidoutProb.write('id Mmin window Probability name\n')

        # Cycle for number of faults
        for i in range(nfault):
            magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
            M = 10 ** (c * magnitude_range + d)
            pdf_mag = norm.pdf(magnitude_range, mag[i], sdmag[i])
            total_moment = np.sum(pdf_mag * M)
            ratio = Morate[i] / total_moment
            balanced_pdf_moment = ratio * pdf_mag
            Mo_balanced = np.sum(balanced_pdf_moment * M)  # for a check
            CHgaussRATES = balanced_pdf_moment
            cumCHgaussRATES = np.flip(np.cumsum(np.flip(CHgaussRATES)))

            Mag_min = magnitude_range[0]
            out_Rates = [id[i], Mag_min, bin] + CHgaussRATES.tolist()
            out_Prob = [id[i], Mag_min, w, Hpois[i]]

            fault_data = {}
            fault_data = {
                "rates": out_Rates # Your 'i' index
            }

            # Adding the outputs of the moment budget to the faults dictionary
            faults[faultname[i]].update(fault_data)


            # Writing to output files
            rates_str = ', '.join(f"{rate:5.4e}" for rate in CHgaussRATES)
            fidout.write(f"{id[i]}, {Mag_min:3.1f}, {bin:3.1f}, {rates_str}, {faultname[i]}\n")
            fidoutProb.write(f"{id[i]}, {Mag_min:3.1f}, {w}, {Hpois[i]:5.3e}, {faultname[i]}\n")

            # Plotting
            plt.figure(i)
            plt.semilogy(magnitude_range, cumCHgaussRATES, 'ok')
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultname[i])
            figname = f"./output_files/{outputfilename}_AR_ChGaussPoisson_rates_{faultname[i]}"
            plt.savefig(figname, format='epsc')


    export_faults_to_xml(faults, outputfilename)
    aaa=111
    return Mo_balanced


# Example of function usage
# output = CHGaussPoiss(c, d, outputfilename, faultname, mag, sdmag, Morate, id, nfault, w, Hpois, bin)




        # Open the XML file for writing


# Usage example with sample data




def CHGaussBPT(faults, c, d, outputfilename, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin):
    outputname = f"{outputfilename}_AR_ChGaussBPT_rates.txt"
    outputnameProbability = f"{outputfilename}_AR_ChGaussBPT_Probability.txt"

    # Create output_files directory if it doesn't exist
    os.makedirs('./output_files/', exist_ok=True)

    # Open a file for writing the output
    with open(f'./output_files/{outputname}', 'w') as fidout, open(f'./output_files/{outputnameProbability}',
                                                                   'w') as fidoutProb:
        # Print a title, followed by a blank line
        fidout.write('id Mmin bin rates name\n')
        fidoutProb.write('id Mmin window Probability name\n')

        # Calculate a fictious Tmean following Pace et al., 2006 and a fictious Mo rate
        Tfict = (-1 * w) / np.log(1 - Hbpt)
        Morate_fict = Morate * (Tmean / Tfict)

        Mo_balanced_fict = []

        # Cycle for number of faults
        for i in range(nfault):
            magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
            M = 10 ** (c * magnitude_range + d)
            pdf_mag = norm.pdf(magnitude_range, mag[i], sdmag[i])
            total_moment = np.sum(pdf_mag * M)
            ratio = Morate_fict[i] / total_moment
            balanced_pdf_moment = ratio * pdf_mag
            Mo_balanced_fict.append(np.sum(balanced_pdf_moment * M))
            CHgaussRATES = balanced_pdf_moment
            cumCHgaussRATES = np.flip(np.cumsum(np.flip(CHgaussRATES)))

            Mag_min = magnitude_range[0]
            out_Rates = [id[i], Mag_min, bin] + CHgaussRATES.tolist()
            out_Prob = [id[i], Mag_min, w, Hbpt[i]]

            fault_data = {}
            # Assuming 'i' is an index variable
            # You can replace this with your actual index
            # Create a dictionary to store fault data
            fault_data = {
                "rates": out_Rates # Your 'i' index
            }

            # Adding the outputs of the moment budget to the faults dictionary
            faults[faultname[i]].update(fault_data)

            # Writing to output files
            fidout.write(f"{id[i]}, {Mag_min:3.1f}, {bin:3.1f}, " + ', '.join(
                f"{rate:5.4e}" for rate in CHgaussRATES) + f", {faultname[i]}\n")
            fidoutProb.write(f"{id[i]}, {Mag_min:3.1f}, {w}, {Hbpt[i]:5.3e}, {faultname[i]}\n")

            # Plotting
            plt.figure(i)
            plt.semilogy(magnitude_range, cumCHgaussRATES, 'ok')
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultname[i])
            figname = f"./output_files/{outputfilename}_AR_ChGaussBPT_rates_{faultname[i]}"
            plt.savefig(figname, format='eps')

    export_faults_to_xml(faults, outputfilename)


    aaaaa=111

    return Mo_balanced_fict


# Example of function usage
# output = CHGaussBPT(c, d, outputfilename, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin)


def kin2coeff(ScR):
    # Check if ScR is a list, and if so, convert it to a string
    if isinstance(ScR, list):
        ScR = ''.join(ScR)

    # Initialize coefficient matrices
    coeff = None
    ARtable = None

    if ScR[:4].upper() == 'WC94':
        # Coefficients by Wells & Coppersmith, 1994
        # Coefficients to calculate Magnitude from length and area
        # Matrix format: [aRLD, bRLD, sdRLD, aRA, bRA, sdRA]
        coeff = [[4.34, 1.54, 0.31, 3.93, 1.02, 0.25],
                 [4.49, 1.49, 0.26, 4.33, 0.90, 0.25],
                 [4.33, 1.49, 0.24, 3.98, 1.02, 0.23],
                 [4.38, 1.49, 0.26, 4.07, 0.98, 0.24]]

    elif ScR[:4].upper() == 'LE10':
        # Coefficients by Leonard, 2010
        # Coefficients to calculate Moment from length and area
        # Matrix format: [aRLD, bRLDmin, bRLDmax, aRA, bRAmin, bRAmax]
        coeff = [[2.5, 7.53, 8.51, 1.5, 5.69, 6.6],
                 [1.5, 12.01, 12.88, 1.5, 5.69, 6.47],
                 [2.5, 7.87, 8.28, 1.5, 6.22, 6.52]]

    elif ScR[:4].upper() == 'AZ15' or ScR[:4].upper() == 'VOLC':
        # Coefficients by D'Amico and Azzaro, 2014 and Villamor 2001
        # Coefficients to calculate Ml (D'Amico) and Mw (Villamor) from
        # surface rupture length (D'Amico) and RA (Villamor)
        # Matrix format: [aSRLmin, aSRLmax, bSRLmin, bSRLmax]
        coeff = [[3.239, 3.543, 1.662, 2.49, 3.39, 1.33]]

    # ASPECT RATIO coefficients by Pace et al., 2002 (BGTA)
    # Matrix format: [aAS, bAS, sdAS]
    # These coefficients are ALWAYS used for computing Magnitudes
    ARtable = [[3.0939, 1.2501, 0.25],
               [-4.4543, 2.1992, 0.25],
               [-7.096, 2.9807, 0.25],
               [-2.3725, 1.9354, 0.25]]

    return coeff, ARtable

################################### Coefficient to magnitude function

def coeff2mag(ScR, coeff, Length, Width, ARtable, mu, straindrop):
    print("Hi")
    print(Length)
    # Compute Magnitudes from scale- and from Aspect Ratio-relationship
    # Check if ScR is a list, and if so, convert it to a string
    if isinstance(ScR, list):
        ScR = ''.join(ScR)
    # Initialize variables
    MRLD = None
    MRA = None
    dMRLD = None
    dMRA = None
    MAR = None
    ar_coeff = None
    LAR = None
    legends_Mw = None

    # Check for the kind of faulting

    if ScR.upper() == 'WC94-N':
        wc = coeff[0]
        ar_coeff = ARtable[0]
        Length_km = Length / 1000
        Width_km = Width / 1000

        MRLD = wc[0] + wc[1] * math.log10(Length_km)
        MRA = wc[3] + wc[4] * math.log10(Length_km * Width_km)
        dMRLD = wc[2]
        dMRA = wc[5]



        legends_Mw = ['MRLD', 'MRA']

    elif ScR.upper() == 'WC94-R':
        wc = coeff[1]
        ar_coeff = ARtable[1]
        Length_km = Length / 1000
        Width_km = Width / 1000

        MRLD = wc[0] + wc[1] * math.log10(Length_km)
        MRA = wc[3] + wc[4] * math.log10(Length_km * Width_km)
        dMRLD = wc[2]
        dMRA = wc[5]

        legends_Mw = ['MRLD', 'MRA']

    elif ScR.upper() == 'WC94-S':
        wc = coeff[2]
        ar_coeff = ARtable[2]
        Length_km = Length / 1000
        Width_km = Width / 1000

        MRLD = wc[0] + wc[1] * math.log10(Length_km)
        MRA = wc[3] + wc[4] * math.log10(Length_km * Width_km)
        dMRLD = wc[2]
        dMRA = wc[5]

        legends_Mw = ['MRLD', 'MRA']

    elif ScR.upper() == 'WC94-A':
        wc = coeff[3]
        ar_coeff = ARtable[3]
        Length_km = Length / 1000
        Width_km = Width / 1000

        MRLD = wc[0] + wc[1] * math.log10(Length_km)
        MRA = wc[3] + wc[4] * math.log10(Length_km * Width_km)
        dMRLD = wc[2]
        dMRA = wc[5]

        legends_Mw = ['MRLD', 'MRA']

    # LEONARD 2010 EQUATIONS
    elif ScR.upper() in ['LE10-N', 'LE10-R', 'LE10-D']:
        leo = coeff[0]
        ar_coeff = ARtable[3]

        MRLDmin = (2 / 3) * math.log10(10**(leo[1] + leo[0] * math.log10(Length))) - 6.07
        MRAmin = (2 / 3) * math.log10(10 ** (leo[4] + leo[3] * math.log10(Length * Width))) - 6.07
        MRLDmax = (2 / 3) * math.log10(10**(leo[2] + leo[0] * math.log10(Length))) - 6.07
        MRAmax = (2 / 3) * math.log10(10**(leo[5] + leo[3] * math.log10(Length * Width))) - 6.07
        MRLD = MRLDmin + ((MRLDmax - MRLDmin) / 2)
        MRA = MRAmin + ((MRAmax - MRAmin) / 2)
        dMRLD = (MRLDmax - MRLDmin) / 2
        dMRA = (MRAmax - MRAmin) / 2

        legends_Mw = ['MRLD', 'MRA']

    elif ScR.upper() == 'LE10-S':
        leo = coeff[1]
        ar_coeff = ARtable[2]

        MRLDmin = (2 / 3) * math.log10(10**(leo[1] + leo[0] * math.log10(Length))) - 6.07
        MRAmin = (2 / 3) * math.log10(10 ** (leo[4] + leo[3] * math.log10(Length * Width))) - 6.07
        MRLDmax = (2 / 3) * math.log10(10**(leo[2] + leo[0] * math.log10(Length))) - 6.07
        MRAmax = (2 / 3) * math.log10(10**(leo[5] + leo[3] * math.log10((Length) * (Width)))) - 6.07
        MRLD = MRLDmin + ((MRLDmax - MRLDmin) / 2)
        MRA = MRAmin + ((MRAmax - MRAmin) / 2)
        dMRLD = (MRLDmax - MRLDmin) / 2
        dMRA = (MRAmax - MRAmin) / 2

        legends_Mw = ['MRLD', 'MRA']

    elif ScR.upper() in ['LE10-SCR', 'LE10-STABLE']:
        leo = coeff[2]
        ar_coeff = ARtable[3]

        MRLDmin = (2 / 3) * math.log10(10**(leo[1] + leo[0] * math.log10(Length))) - 6.07
        MRAmin = (2 / 3) * math.log10(10**(leo[4] + leo[3] * math.log10(Length * Width)) - 6.07)
        MRLDmax = (2 / 3) * math.log10(10**(leo[2] + leo[0] * math.log10(Length))) - 6.07
        MRAmax = (2 / 3) * math.log10(10**(leo[5] + leo[3] * math.log10(Length * Width)) - 6.07)
        MRLD = MRLDmin + ((MRLDmax - MRLDmin) / 2)
        MRA = MRAmin + ((MRAmax - MRAmin) / 2)
        dMRLD = (MRLDmax - MRLDmin) / 2
        dMRA = (MRAmax - MRAmin) / 2

        legends_Mw = ['MRLD', 'MRA']

    # AZZARO et al. 2014 and VILLAMOR 2001 EQUATIONS (VOLCANIC CONTEXT NORMAL KIN)
    elif ScR.upper() == 'VOLC' or ScR[:4].upper() == 'AZ15':
        vol = coeff[0]
        ar_coeff = ARtable[0]
        Length_km = Length / 1000
        Width_km = Width / 1000

        Mlmin = vol[0] + vol[2] * math.log10(Length_km)
        Mlmax = vol[1] + vol[3] * math.log10(Length_km)
        Mwmin = (vol[4] + vol[5] * math.log10(Length_km * Width_km)) - 0.195
        Mwmax = (vol[4] + vol[5] * math.log10(Length_km * Width_km)) + 0.195
        MRLD = Mlmin + ((Mlmax - Mlmin) / 2)  # Ml
        MRA = Mwmin + ((Mwmax - Mwmin) / 2)  # Mw
        dMRLD = (Mlmax - Mlmin) / 2  # dMl
        dMRA = (Mwmax - Mwmin) / 2  # dMw

        legends_Mw = ['MlDA', 'MwVi']

    # Aspect Ratio Control Formula (Pace and Peruzza, 2002)
    # Note that here Length and Width are expressed in km

    Width_km = Width / 1000

    LAR = ar_coeff[0] + ar_coeff[1] * (Width_km)

    # Check if Length from Aspect Ratio is not greater than input length
    # Now dimensions are in meters


    LAR = LAR * 1000

    # Calculate the moment magnitude from LAR
    MAR = (2 / 3) * (math.log10(straindrop * mu * LAR ** 2 * Width) - 9.05)

    return MRLD, MRA, dMRLD, dMRA, MAR, ar_coeff, LAR, legends_Mw




# Conflation function


def conflate_pdfs(x, pdfs):
    # Initialize the conflated distribution as a uniform distribution
    conflated = np.ones(x.shape) / (x[-1] - x[0])

    for i in range(pdfs.shape[0]):
        pdf = pdfs[i, :] / np.trapz(x, pdfs[i, :])
        conflated = conflated * pdf

    # Normalize the conflated distribution
    conflated = -conflated / np.trapz(x, np.abs(conflated))
    conflated=np.abs(conflated) # Ensure positive values

    return conflated

# # Functions to apply MFD options
# # Characteristic Gaussian (BPT)
#
# def CHGaussBPT(c, d, outputfilename, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin):
#     outputname = f'{outputfilename}_AR_ChGaussBPT_rates.txt'
#     outputnameProbability = f'{outputfilename}_AR_ChGaussBPT_Probability.txt'
#
#     # Open files for writing the output
#     with open(f'./output_files/{outputname}', 'w') as fidout, open(f'./output_files/{outputnameProbability}', 'w') as fidoutProb:
#         # Write headers
#         fidout.write('id Mmin bin rates name\n')
#         fidoutProb.write('id Mmin window Probability name\n')
#
#         # Calculate fictitious Tmean following Pace et al., 2006 and a fictitious Mo rate
#         Tfict = (-1 * w) / np.log(1 - Hbpt)
#         Morate_fict = Morate * (Tmean / Tfict)
#
#         for i in range(nfault):
#             magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
#             M = 10 ** (c * magnitude_range + d)
#
#             pdf_mag = stats.norm.pdf(magnitude_range, loc=mag[i], scale=sdmag[i])
#             total_moment = np.sum(pdf_mag * M)
#
#             ratio = Morate_fict[i] / total_moment
#             balanced_pdf_moment = ratio * pdf_mag
#
#             Mo_balanced_fict = np.sum(balanced_pdf_moment * M)
#
#             CHgaussRATES = balanced_pdf_moment
#             cumCHgaussRATES = np.cumsum(np.flipud(CHgaussRATES))
#             Mag_min = magnitude_range[0]
#
#             out_Rates = [id[i], Mag_min, bin] + CHgaussRATES.tolist()
#             out_Prob = [id[i], Mag_min, w, Hbpt[i]]
#
#             fidout.write(f"{out_Rates[0]}, {out_Rates[1]:.1f}, {out_Rates[2]:.1f},")
#             fidout.write(" ")
#             fidout.write(" ".join([f"{rate:.4e}" for rate in out_Rates[3:]]))
#             fidout.write(" ")
#             fidout.write(f",{faultname[i]}\n")
#
#             fidoutProb.write(f"{out_Prob[0]}, {out_Prob[1]:.1f}, {out_Prob[2]}, {out_Prob[3]:.3e},")
#             fidoutProb.write(" ")
#             fidoutProb.write(f"{faultname[i]}\n")
#
#             # Create a figure for each fault
#             plt.figure(i)
#             plt.semilogy(magnitude_range, cumCHgaussRATES, 'ok')
#             fault = faultname[i]
#             figname = f'./output_files/{outputfilename}_AR_ChGaussBPT_rates_{fault}.png'
#             plt.xlabel('magnitude')
#             plt.ylabel('annual cumulative rates')
#             plt.title(fault)
#             plt.savefig(figname)
#
#
# # Characteristic Gaussian (Poiss)
#
# def CHGaussPoiss(c, d, outputfilename, faultname, mag, sdmag, Morate, id, nfault, w, Hpois, bin):
#     outputname = f"{outputfilename}_AR_ChGaussPoisson_rates.txt"
#     outputnameProbability = f"{outputfilename}_AR_ChGaussPoisson_Probability.txt"
#
#     # Open two files for writing the outputs
#     with open(f'./output_files/{outputname}', 'w') as fidout, open(f'./output_files/{outputnameProbability}',
#                                                                    'w') as fidoutProb:
#         # Print a title, followed by a blank line
#         fidout.write('id Mmin bin rates name\n')
#         fidoutProb.write('id Mmin window Probability name\n')
#
#         for i in range(nfault):
#             magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
#             M = 10 ** (c * magnitude_range + d)
#
#             pdf_mag = stats.norm.pdf(magnitude_range, loc=mag[i], scale=sdmag[i])
#             total_moment = np.sum(pdf_mag * M)
#             ratio = Morate[i] / total_moment
#             balanced_pdf_moment = ratio * pdf_mag
#             Mo_balanced = np.sum(balanced_pdf_moment * M)
#
#             CHgaussRATES = balanced_pdf_moment
#             cumCHgaussRATES = np.cumsum(CHgaussRATES)
#             Mag_min = magnitude_range[0]
#             out_Rates = [id[i], Mag_min, bin] + list(CHgaussRATES)
#             out_Prob = [id[i], Mag_min, w, Hpois[i]]
#
#             # Write output to files
#             fidout.write(f'{id[i]}, {Mag_min:.1f}, {bin:.1f}, ')
#             fidout.write(' '.join([f'{rate:.4e}' for rate in CHgaussRATES]) + ', ')
#             fidout.write(f'{faultname[i]}\n')
#
#             fidoutProb.write(f'{id[i]}, {Mag_min:.1f}, {w}, {Hpois[i]:.3e}, {faultname[i]}\n')
#
#             # Plot figures and save output files
#             plt.figure(i)
#             plt.semilogy(magnitude_range, cumCHgaussRATES, 'ok')
#             plt.xlabel('magnitude')
#             plt.ylabel('annual cumulative rates')
#             plt.title(faultname[i])
#             figname = f'./output_files/{outputfilename}_AR_ChGaussPoisson_rates_{faultname[i]}'
#             plt.savefig(figname + '.png')
#             plt.close(i)
#
#
#
# # Characteristic Gaussian (Prob)
#
# def CHGaussProb(c, d, outputfilename, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, ProbabilityOfOccurrence,
#                 bin):
#     outputname = f"{outputfilename}_AR_ChGaussProb.txt"
#
#     # Open a file for writing the output
#     with open(f'./output_files/{outputname}', 'w') as fidout:
#         # Print a title, followed by a blank line
#         fidout.write('id Mmin bin rates name\n')
#
#         # Calculate a fictitious Tmean following Pace et al., 2006
#         Tfict = (-1 * w) / np.log(1 - ProbabilityOfOccurrence)
#         Morate_fict = Morate * (Tmean / Tfict)
#
#         for i in range(nfault):
#             magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
#             M = 10 ** (c * magnitude_range + d)
#
#             pdf_mag = stats.norm.pdf(magnitude_range, loc=mag[i], scale=sdmag[i])
#             total_moment = np.sum(pdf_mag * M)
#
#             ratio = Morate_fict[i] / total_moment
#             balanced_pdf_moment = ratio * pdf_mag
#
#             Mo_balanced_fict = np.sum(balanced_pdf_moment * M)
#
#             CHgaussRATES = balanced_pdf_moment
#             cumCHgaussRATES = np.cumsum(CHgaussRATES)
#             Mag_min = magnitude_range[0]
#             out_Rates = [id[i], Mag_min, bin] + list(CHgaussRATES)
#
#             # Write output to the file
#             fidout.write(f'{id[i]}, {Mag_min:.1f}, {bin:.1f}, ')
#             fidout.write(' '.join([f'{rate:.4e}' for rate in CHgaussRATES]) + ', ')
#             fidout.write(f'{faultname[i]}\n')
#
#             plt.figure(i)
#             plt.semilogy(magnitude_range, cumCHgaussRATES, 'ok')
#             fault = faultname[i]
#             figname = f'./output_files/{outputfilename}_AR_ChGaussProb_rates_{fault}'
#             plt.xlabel('magnitude')
#             plt.ylabel('annual cumulative rates')
#             plt.title(fault)
#             plt.savefig(figname + '.png')
#             plt.close(i)
#
#
# # Truncated Gutenberg Richter
#
# def TruncatedGR(c, d, outputfilename, faultname, mag, mt, Morate, id, nfault, bin, b):
#     outputname = f"{outputfilename}_AR_TruncatedGR.txt"
#
#     with open(f"./output_files/{outputname}", "w") as fidout:
#         fidout.write("id Mmin bin rates name\n")
#
#         DATA = {}
#
#         for i in range(nfault):
#             magnitude_range = np.arange(mt[i], mag[i] + bin, bin)
#             M = 10 ** (c * magnitude_range + d)
#             Beta = (2 / 3) * b[i]
#             Mt = 10 ** (c * mt[i] + d)
#             Mxp = 10 ** (c * (mag[i] + bin) + d)
#             TruncGR = (((Mt / M) ** Beta - (Mt / Mxp) ** Beta) / (1 - (Mt / Mxp) ** Beta))
#             Incremental = np.concatenate((np.diff(TruncGR[::-1])[::-1], [TruncGR[-1]]))
#             Incremental_Morate = Incremental * M
#             Incremental_Morate_balanced = (Incremental_Morate * Morate[i]) / np.sum(Incremental_Morate)
#             cons_tassi_ind = (Incremental * Incremental_Morate_balanced) / Incremental_Morate
#             cumulative_rates = np.cumsum(cons_tassi_ind[::-1])[::-1]
#             Mbalanced = np.sum(cons_tassi_ind * M)
#             out_Rates = [id[i], mt[i], bin, cons_tassi_ind]
#             DATA[faultname[i]] = {
#                 "rates": cons_tassi_ind.tolist(),
#                 "bin": bin
#             }
#             fidout.write(f"{out_Rates[0]:d} {out_Rates[1]:.1f} {out_Rates[2]:.1f} ")
#             fidout.write(" ".join([f"{rate:.4e}" for rate in out_Rates[3]]))
#             fidout.write(f" {faultname[i]}\n")
#
#             plt.figure(i)
#             plt.semilogy(magnitude_range, cumulative_rates, 'ok')
#             fault = faultname[i]
#             figname = f"./output_files/{outputfilename}_AR_TruncatedGR_rates_{fault}.png"
#             plt.xlabel('magnitude')
#             plt.ylabel('annual cumulative rates')
#             plt.title(fault)
#             plt.savefig(figname)
#
#         export_faults_to_xml(DATA)
#
#
# def export_faults_to_xml(DATA):
#     root = ET.Element("Faults")
#
#     for faultname, data in DATA.items():
#         fault_elem = ET.SubElement(root, "Fault")
#         fault_elem.set("name", faultname)
#
#         rates_elem = ET.SubElement(fault_elem, "Rates")
#         for rate in data["rates"]:
#             rate_elem = ET.SubElement(rates_elem, "Rate")
#             rate_elem.text = str(rate)
#
#         bin_elem = ET.SubElement(fault_elem, "Bin")
#         bin_elem.text = str(data["bin"])
#
#     tree = ET.ElementTree(root)
#     tree.write("faults.xml")
#
#
# if __name__ == "__main__":
#     outputfilename = "output"
#     mag = np.array([5.0, 6.0, 7.0])
#     id = np.array([1, 2, 3])
#     nfault = len(mag)
#     faultname = ["Fault1", "Fault2", "Fault3"]
#     mt = np.array([5.5, 6.5, 7.5])
#     Morate = np.array([0.1, 0.2, 0.3])
#     bin = 0.1
#     b = np.array([0.1, 0.2, 0.3])
#
#     TruncatedGR(c, d, outputfilename, faultname, mag, mt, Morate, id, nfault, bin, b)
#
#
#
#
# # truncGaussDistObsMag
#
# def truncGaussDistObsMag(pdf_magnitudes, x_range_of_mag, mag, sdmag):
#     dist_mag = pdf_magnitudes[-1]
#     dist_obs_mag = np.zeros(dist_mag.size)
#     lower_threshold = mag - sdmag
#     upper_threshold = mag + sdmag
#     lowest_value = np.argmax(x_range_of_mag <= lower_threshold)
#     highest_value = np.argmax(x_range_of_mag <= upper_threshold)
#
#     dist_obs_mag[lowest_value:highest_value+1] = dist_mag[lowest_value:highest_value+1]
#
#     return dist_obs_mag
#
# ##########################################################################################
#
# # # Logfile Ceation for MB
# #
# #
# # def LogFileMB(inputdata, outputfile, data, ScR, faultname, shearmodulusfrominput, straindropfrominput,
# #               M_dM_Lengths_forLOGoutput, out, n_sigma_trunc_magnitudes, yfc):
# #     if os.name == 'nt':
# #         inputdata = inputdata.replace('/', '\\')
# #
# #     with open('./output_files/LogFileMB.txt', 'a') as fidout:
# #         fidout.write('*************************\n\n\n')
# #         appending_pdf = open('./output_files/temporary_pdf_storage.txt', 'r')
# #
# #         # Print the input file name followed by a blank line
# #         fidout.write(f'input file used: {inputdata}\n')
# #         # Print the output file name followed by a blank line
# #         fidout.write(f'output file created: {outputfile}\n')
# #         # Print the time at which you did it followed by a blank line
# #         t = time.localtime()
# #         fidout.write('executed at: ')
# #         fidout.write(f'{t.tm_year} {t.tm_mon} {t.tm_mday} {t.tm_hour} {t.tm_min} {t.tm_sec}\n\n')
# #
# #         # Print the input data followed by a blank line
# #         fidout.write('input data used for calculation:\n')
# #         # Print a title, followed by a blank line
# #         fidout.write(
# #             'name ScaleRel length Dip Seismogenic_Thickness slipratemin slipratemax Mobs sdMobs occurrenceTime ShearMod(*10^10) StrainDrop(*10^-5) \n')
# #         for i in range(data.shape[0]):
# #             fidout.write(f'{faultname[i]} ')
# #             fidout.write(f'{ScR[i]} ')
# #             fidout.write(
# #                 f'{data[i, 2]:.2f} {data[i, 3]:.2f} {data[i, 4]:.2f} {data[i, 5]:.2f} {data[i, 6]:.2f} {data[i, 7]:.1f} {data[i, 8]:.1f} ')
# #             fidout.write(f'{data[i, 9]} {shearmodulusfrominput} {straindropfrominput:.5e}\n')
# #         fidout.write('\n\n')
# #
# #         # Print the option used for calculation
# #         fidout.write(f'number of sigma for trucated magnitudes function: {n_sigma_trunc_magnitudes}\n')
# #         fidout.write(f'year for calculating elasped time: {yfc}\n\n')
# #
# #         # Print the output data followed by a blank line
# #         fidout.write('output data obtained:\n')
# #         # Print a title, followed by a blank line
# #         fidout.write(
# #             'name input-Lenght Lenght-from-AspectRatio MMo Mar Mscr1 Mscr2 Mobs dMMo dMar dMscr1 dMscr2 dMobs id Mmax sdMmax Tmean CV Telap Mo-rate\n')
# #         for i in range(out.shape[0]):
# #             fidout.write(f'{faultname[i]} ')
# #             fidout.write(
# #                 f'{out[i, 0]:.2f} {out[i, 1]:.2f} {out[i, 2]:.2f} {out[i, 3]:.2f} {out[i, 4]:.2f} {out[i, 5]:.2f} ')
# #             fidout.write(
# #                 f'{out[i, 6]:.2f} {out[i, 7]:.2f} {out[i, 8]:.2f} {out[i, 9]:.2f} {out[i, 10]:.2f} {out[i, 11]:.2f} ')
# #             fidout.write(
# #                 f'{out[i, 12]} {out[i, 13]:.1f} {out[i, 14]:.1f} {out[i, 15]} {out[i, 16]:.1f} {out[i, 17]} {out[i, 18]:.5e}\n')
# #
# #         # Weights (NaN if unweighted), x_range of magnitude, and their unweighted-values
# #         fidout.write('\n')
# #         fidout.write('weigths(NaN if unweighted, 0 if not used), x_range of magnitude and pdf_of_mag-values:\n\n')
# #         for tline in appending_pdf:
# #             fidout.write(tline)
# #
# #         fidout.write('*************************\n\n\n')
# #
# #
# # LogFileMB(inputdata, outputfile, data, ScR, faultname, shearmodulusfrominput, straindropfrominput,
# #           M_dM_Lengths_forLOGoutput, out, n_sigma_trunc_magnitudes, yfc)
# #
# # # Logfile creation of AR
# #
# # def LogFileAR(inputdata, outputname, data, faultname, Fault_behaviour, MFD, window, binstep, grinputoptions,
# #               ParametersGR):
# #     if os.name == 'nt':
# #         inputdata = inputdata.replace('/', '\\')
# #
# #     with open('./output_files/LogFileAR.txt', 'a') as fidout:
# #         fidout.write('*************************\n\n\n')
# #
# #         # Print the input file name followed by a blank line
# #         fidout.write(f'input file used: {inputdata}\n')
# #         # Print the output file name followed by a blank line
# #         fidout.write(f'output file created: {outputname}\n')
# #         # Print the time at which you did it followed by a blank line
# #         t = time.localtime()
# #         fidout.write('executed at: ')
# #         fidout.write(f'{t.tm_year} {t.tm_mon} {t.tm_mday} {t.tm_hour} {t.tm_min} {t.tm_sec:.2f}\n\n')
# #
# #         # Print the input data followed by a blank line
# #         fidout.write('input data used for calculation:\n')
# #         # Print a title, followed by a blank line
# #         if Fault_behaviour in [1, 2, 4, 5, 7, 8]:
# #             fidout.write('id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) name\n')
# #             for i in range(data.shape[0]):
# #                 fidout.write(
# #                     f'{data[i, 0]} {data[i, 1]:.2f} {data[i, 2]:.2f} {data[i, 3]} {data[i, 4]:.2f} {data[i, 5]} {data[i, 6]:.3e} ')
# #                 fidout.write(f'{faultname[i]} {faultname[i]}')
# #                 fidout.write('\n')
# #         elif Fault_behaviour in [3, 6]:
# #             fidout.write('id Mchar sdMchar Tmean alfa Telap Mo_rate(N/m2) Probability name\n')
# #             for i in range(data.shape[0]):
# #                 fidout.write(
# #                     f'{data[i, 0]} {data[i, 1]:.2f} {data[i, 2]:.2f} {data[i, 3]} {data[i, 4]:.2f} {data[i, 5]} {data[i, 6]:.3e} {data[i, 7]:.3e} ')
# #                 fidout.write(f'{faultname[i]} {faultname[i]}')
# #                 fidout.write('\n')
# #
# #         fidout.write('\n\n')
# #
# #         # Print the options used for calculation
# #         fidout.write(f'Magnitude-Frequency-Distribution: {MFD}\n')
# #         if not math.isnan(grinputoptions):
# #             fidout.write('GR options (id, b-value, Mt):\n')
# #             for i in range(len(ParametersGR)):
# #                 fidout.write(f'{ParametersGR[i, 0]} {ParametersGR[i, 1]:.2f} {ParametersGR[i, 2]:.2f}\n')
# #         fidout.write(f'forecast period: {window}\n')
# #         fidout.write(f'binstep: {binstep}\n\n')
# #
# #         fidout.write('*************************\n\n\n')
# #
# #
# # # Example usage
# # inputdata = "C:/path/to/input/file"
# # outputname = "LogFileAR_output.txt"
# # data = [[1, 2.5, 0.1, 50, 1.2, 100, 1.5], [2, 3.0, 0.2, 45, 1.3, 150, 1.6]]
# # faultname = ["Fault1", "Fault2"]
# # Fault_behaviour = 1
# # MFD = "SomeMFD"
# # window = 10
# # binstep = 1
# # grinputoptions = 0
# # ParametersGR = [[1, 2.1, 4.5], [2, 2.2, 4.6]]
# # LogFileAR(inputdata2, outputname, data, faultname, Fault_behaviour, MFD, window, binstep, grinputoptions, ParametersGR)
#
# ##########################################################################################
#
#
# # kin2coeff
#
# def kin2coeff(ScR):
#     if ScR.lower().startswith('wc94'):
#         # Coefficients by Wells & Coppersmith, 1994
#         # Coefficients to calculate Magnitude from length and area
#         # Matrix: row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for ALL
#         # Each row contains aRLD, bRLD, sdRLD, aRA, bRA, sdRA
#         coeff = [
#             [4.34, 1.54, 0.31, 3.93, 1.02, 0.25],
#             [4.49, 1.49, 0.26, 4.33, 0.90, 0.25],
#             [4.33, 1.49, 0.24, 3.98, 1.02, 0.23],
#             [4.38, 1.49, 0.26, 4.07, 0.98, 0.24]
#         ]
#     elif ScR.lower().startswith('le10'):
#         # Coefficients by Leonard, 2010
#         # Coefficients to calculate Moment from length and area
#         # Matrix: row 1 for DS dip slip; row 2 for SS strike slip; row 3 for SCR Stable Continental Regions
#         # Each row contains aRLD, bRLDmin, bRLDmax, aRA, bRAmin, bRAmax
#         # RLD=L RA=A (see Table 5 in Leonard, 2010)
#         # DS: A>5,500m2 SS:A45,000 m2 SCR: A>2,500 m2
#         coeff = [
#             [2.5, 7.53, 8.51, 1.5, 5.69, 6.6],
#             [1.5, 12.01, 12.88, 1.5, 5.69, 6.47],
#             [2.5, 7.87, 8.28, 1.5, 6.22, 6.52]
#         ]
#     elif ScR.lower().startswith('az15') or ScR.lower().startswith('volc'):
#         # Coefficients by D'Amico and Azzaro, 2014 and Villamor 2001
#         # Coefficients to calculate Ml (D'Amico) and Mw (Villamor) from surface rupture length (D'Amico) and RA (Villamor)
#         # Matrix: row 1 for D'Amico; row 2 for Villamor
#         # Each row contains aSRLmin, aSRLmax, bSRLmin, bSRLmax, SRL=surface rupture length
#         coeff = [
#             [3.239, 3.543, 1.662, 2.49],
#             [3.39, 1.33]
#         ]
#
#     # Aspect Ratio coefficients by Pace et al., 2002 (BGTA)
#     # Matrix: row 1 for normal; row 2 for reverse; row 3 for strikeslip; row 4 for ALL DIP
#     # Each row contains aAS, bAS, sdAS
#     # These coefficients are ALWAYS used for computing Magnitudes
#     ARtable = [
#         [3.0939, 1.2501, 0.25],
#         [-4.4543, 2.1992, 0.25],
#         [-7.096, 2.9807, 0.25],
#         [-2.3725, 1.9354, 0.25]
#     ]
#
#     return coeff, ARtable
#
#
# # Example usage:
# ScR = 'WC94'
# coeff, ARtable = kin2coeff(ScR)
# print("Coefficients:")
# for row in coeff:
#     print(row)
#
# print("\nAspect Ratio Coefficients:")
# for row in ARtable:
#     print(row)
#
#
# ########## coeff2mag
# def coeff2mag(ScR, coeff, Length, Width, ARtable, mu, straindrop):
#     # Coefficients for computing magnitude from length and area
#     coeff = np.array(coeff)
#     Length_km = Length / 1000
#     Width_km = Width / 1000
#
#     if ScR.lower() in ['wc94-n', 'wc94-r', 'wc94-s', 'wc94-a']:
#         wc = coeff[ScR.lower().endswith(('-n', '-r', '-s', '-a'))]
#         ar_coeff = ARtable[ScR.lower().endswith(('-n', '-r', '-s', '-a'))]
#
#         MRLD = wc[0] + wc[1] * np.math.log10(Length_km)
#         MRA = wc[3] + wc[4] * np.math.log10(Length_km * Width_km)
#         dMRLD = wc[2]
#         dMRA = wc[5]
#         legends_Mw = ['MRLD', 'MRA']
#
#     elif ScR.lower() in ['le10-n', 'le10-r', 'le10-d', 'le10-s', 'le10-scr', 'le10-stable']:
#         leo = coeff[ScR.lower().endswith(('-n', '-r', '-d', '-s', '-scr', '-stable'))]
#         ar_coeff = ARtable[3]
#
#         MRLDmin = (2 / 3) * np.math.log10(10 ** (leo[1] + leo[0] * np.math.log10(Length))) - 6.07
#         MRAmin = (2 / 3) * np.math.log10(10 ** (leo[4] + leo[3] * np.math.log10((Length) * (Width))) - 6.07)
#         MRLDmax = (2 / 3) * np.math.log10(10 ** (leo[2] + leo[0] * np.math.log10(Length))) - 6.07
#         MRAmax = (2 / 3) * np.math.log10(10 ** (leo[5] + leo[3] * np.math.log10((Length) * (Width))) - 6.07)
#         MRLD = MRLDmin + ((MRLDmax - MRLDmin) / 2)
#         MRA = MRAmin + ((MRAmax - MRAmin) / 2)
#         dMRLD = (MRLDmax - MRLDmin) / 2
#         dMRA = (MRAmax - MRAmin) / 2
#         legends_Mw = ['MRLD', 'MRA']
#
#     elif ScR.lower() in ['az15', 'volc']:
#         vol = coeff[ScR.lower().startswith(('az15', 'volc'))]
#         ar_coeff = ARtable[0]
#         Length_km = Length / 1000
#         Width_km = Width / 1000
#
#         Mlmin = vol[0] + vol[2] * np.math.log10(Length_km)
#         Mlmax = vol[1] + vol[3] * np.math.log10(Length_km)
#         Mwmin = (vol[4] + vol[5] * np.math.log10(Length_km * Width_km)) - 0.195
#         Mwmax = (vol[4] + vol[5] * np.math.log10(Length_km * Width_km)) + 0.195
#         MRLD = Mlmin + ((Mlmax - Mlmin) / 2)  # Ml
#         MRA = Mwmin + ((Mwmax - Mwmin) / 2)  # Mw
#         dMRLD = (Mlmax - Mlmin) / 2  # dMl
#         dMRA = (Mwmax - Mwmin) / 2  # dMw
#         legends_Mw = ['MlDA', 'MwVi']
#
#     # Aspect Ratio Control Formula (Pace and Peruzza, 2002)
#     LAR = ar_coeff[0] + ar_coeff[1] * (Width_km)
#     LAR = min(Length, LAR * 1000)  # Check if Length from Aspect Ratio is not greater than the input length
#     LAR = LAR * 1000  # Now dimensions are in metres
#     # Calculate the moment magnitude from LAR
#     MAR = (2 / 3) * (np.math.log10(straindrop * mu * LAR ** 2 * Width) - 9.05)
#
#     return MRLD, MRA, dMRLD, dMRA, MAR, ar_coeff, LAR, legends_Mw
#
# # Example usage:
# ScR = 'WC94-N'
# Length = 10  # km
# Width = 5  # km
# mu = 3.0
# straindrop = 0.02
# coeff, ARtable = kin2coeff(ScR)
# MRLD, MRA, dMRLD, dMRA, MAR, ar_coeff, LAR, legends_Mw = coeff2mag(ScR, coeff, Length, Width, ARtable, mu, straindrop)
# print("MRLD:", MRLD)
# print("MRA:", MRA)
# print("dMRLD:", dMRLD)
# print("dMRA:", dMRA)
# print("MAR:", MAR)
# print("ar_coeff:", ar_coeff)
# print("LAR:", LAR)
# print("legends_Mw:", legends_Mw)
