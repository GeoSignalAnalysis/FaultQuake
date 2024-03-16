# FaultQuake_functions.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 2023

@authors: Nasrin Tavakolizadeh, Hamzeh Mohammadigheymasi
"""

"""
This module contains funcitons utilized in the FaultQuake workflow
"""


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
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize


def export_faults_to_xml(faults, Project_foldername):
    # Check if './output_files/Fault_OQ' folder exists, if not, create it
    os.makedirs(f'./{Project_foldername}/Fault_OQ', exist_ok=True)

    for fault_name in faults:
        # ScR = faults[fault_name]['ScR']
        # Length = np.append(Length, fault)
        xml_filename = f'./{Project_foldername}/Fault_OQ/{fault_name}.xml'
        current_directory = os.getcwd()
        print(f'Current directory: {current_directory}')

        # Determine XML filename from fault name and save in './Sources' folder
        xml_filename = f'./{Project_foldername}/Fault_OQ/{fault_name}.xml'
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


def TruncatedGR(faults, c, d, Project_foldername, faultnames, mags, mts, Morates, ids, nfault, bin, bs):
    outputname = f"{Project_foldername}_SAR_TruncatedGR.txt"
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
            out_rates = cons_tassi_ind.tolist()

            # cumulative_rates = np.cumsum(np.flipud(cons_tassi_ind))
            # cumulative_rates = np.flipud(cumulative_rates)



            fault_data = {}
            fault_data = {
                "rates": out_rates # Your 'i' index
            }

            # Adding the outputs of the moment budget to the faults dictionary
            faults[faultnames[i]].update(fault_data)

            # Update DATA dictionary
            # DATA[faultnames[i]] = {'id': ids[i], 'rates': cons_tassi_ind, 'bin': bin}

            # Write to output file
            fidout.write(f"{ids[i]}, {mts[i]:3.1f}, {bin:3.1f}, ")
            fidout.write(' '.join(f"{rate:5.4e}" for rate in cons_tassi_ind))

            # Plotting
            plt.figure()
            plt.style.use('dark_background')
            plt.figure(figsize=(8, 6))
            plt.semilogy(magnitude_range, cumulative_rates, 'oy')
            plt.fill_between(magnitude_range, 0, cumulative_rates, color='gold', alpha=0.3)
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultnames[i])

            # Saving the figure
            figname = f"{Project_foldername}_SAR_TruncatedGR_rates_{faultnames[i]}.png"
            plt.savefig(os.path.join('./output_files/Figures', figname), format='png')

            plt.show()  # Show the plot


    export_faults_to_xml(faults, Project_foldername)


    return outputname




def CHGaussPoiss(faults, c, d, Project_foldername, faultname, mag, sdmag, Morate, id, nfault, w, Hpois, bin):
    outputname = f"{Project_foldername}_SAR_ChGaussPoisson_rates.txt"
    outputnameProbability = f"{Project_foldername}_SAR_ChGaussPoisson_Probability.txt"

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
            plt.figure()
            plt.style.use('dark_background')
            plt.figure(figsize=(8, 6))
            plt.semilogy(magnitude_range, cumCHgaussRATES, 'oy')
            plt.fill_between(magnitude_range, 0, cumCHgaussRATES, color='blue', alpha=0.3)
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultname[i])

            # Saving the figure
            figname = f"{Project_foldername}_SAR_TruncatedGR_rates_{faultname[i]}.png"
            plt.savefig(os.path.join('./output_files/Figures', figname), format='png')

            plt.show()  # Show the plot


    export_faults_to_xml(faults, Project_foldername)
    return Mo_balanced


# Example of function usage
# output = CHGaussPoiss(c, d, Project_foldername, faultname, mag, sdmag, Morate, id, nfault, w, Hpois, bin)




        # Open the XML file for writing


# Usage example with sample data




def CHGaussBPT(faults, c, d, Project_foldername, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin):
    outputname = f"{Project_foldername}_SAR_ChGaussBPT_rates.txt"
    outputnameProbability = f"{Project_foldername}_SAR_ChGaussBPT_Probability.txt"

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
            plt.figure()
            plt.style.use('dark_background')
            plt.figure(figsize=(8, 6))
            plt.semilogy(magnitude_range, cumCHgaussRATES, 'oy')
            plt.fill_between(magnitude_range, 0, cumCHgaussRATES, color='magenta', alpha=0.3)
            plt.xlabel('magnitude')
            plt.ylabel('annual cumulative rates')
            plt.title(faultname[i])

            # Saving the figure
            figname = f"{Project_foldername}_SAR_TruncatedGR_rates_{faultname[i]}.png"
            plt.savefig(os.path.join('./output_files/Figures', figname), format='png')

            plt.show()  # Show the plot

            www=353635




    export_faults_to_xml(faults, Project_foldername)

    return Mo_balanced_fict


# Example of function usage
# output = CHGaussBPT(c, d, Project_foldername, faultname, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin)


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

