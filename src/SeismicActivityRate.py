"""
Created on September 2023

@authors: Nasrin Tavakolizadeh, Hamzeh Mohammadigheymasi
"""

"""
Runs Moment of Budget and Activity rate functions
"""

 import math
import numpy
import numpy as np
import os
from scipy.stats import norm
from src.FaultQuake_functions import kin2coeff, coeff2mag, conflate_pdfs, CHGaussBPT, CHGaussPoiss, TruncatedGR
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import invgauss
import scipy.io


import statsmodels.api as sm



def momentbudget(faults, Zeta, Khi, Siggma, ProjFol, logical_nan,  logical_nan_sdmag):
    for fault, values in faults.items():
        if 'ShearModulus' in values and (values['ShearModulus'] is None):
            print(f"Fault {fault} has a 'NaN' ShearModulus value.")
            values['ShearModulus'] = 3 * 1e10
        else:
            values['ShearModulus'] = values['ShearModulus'] * 1e10

        if 'StrainDrop' in values and (values['StrainDrop'] is None):
            values['StrainDrop'] = 3 * 1e-5
        else:
            values['StrainDrop'] = values['StrainDrop'] * 1e-5

        if 'Last_eq_time' in values:
            if values['Last_eq_time'] is None or str(values['Last_eq_time']).lower() == 'nan' or values['Last_eq_time'] == '':
                if logical_nan:
                    values['Last_eq_time'] = math.nan
            else:
                values['Last_eq_time'] = float(values['Last_eq_time'])

        if 'sdmag' in values:
            if values['sdmag'] is None or str(values['sdmag']).lower() == 'nan' or values['sdmag'] == '':
                if logical_nan_sdmag:
                    values['sdmag'] = math.nan
            else:
                values['sdmag'] = float(values['sdmag'])



    kk = 1
    counter = 0
    for fault_name in faults:
        flag_mobs=0
        ScR = faults[fault_name]['ScR']; Length = faults[fault_name]['Length']; Dip = faults[fault_name]['Dip']
        Seismogenic_thickness = faults[fault_name]['Seismogenic_Thickness']; Slipmin = faults[fault_name]['SRmin']
        Slipmax = faults[fault_name]['SRmax']; mag = faults[fault_name]['Mobs']; sdmag = faults[fault_name]['sdMobs']
        Last_eq_time = faults[fault_name]['Last_eq_time']; SCC = faults[fault_name]['SCC']
        ShearModulusFromInputFile = faults[fault_name]['ShearModulus']
        StrainDropFromInputFile = faults[fault_name]['StrainDrop']
        yfc = faults[fault_name]['year_for_calculations']

        if not math.isnan(Last_eq_time):
            Telap = yfc - Last_eq_time  # Perform element-wise subtraction only if Last_eq_time is not NaN


        Telap = yfc - Last_eq_time  # Perform element-wise subtraction
        mu = ShearModulusFromInputFile;         straindrop = StrainDropFromInputFile;         Dip_radians = math.radians(Dip)  # Convert degrees to radians
        sine_Dip = math.sin(Dip_radians);        Length = Length * 1000;        Width = (Seismogenic_thickness * 1000) / sine_Dip
        V = (Slipmin + Slipmax) / 2000;        dV = V - (Slipmin / 1000)
        # DEFINITION OF COEFFICIENTS coefficients by Hanks & Kanamori, 1985
        d = 9.1;   c = 1.5;   dMMO = Siggma  # standard deviation of the magnitude calculated by the seismic moment definition
        # Calculation of the coefficient of the imported scale-relationship
        coeff, ARtable = kin2coeff(ScR)
        MRLD, MRA, dMRLD, dMRA, MAR, ar_coeff, LAR, legends_Mw = coeff2mag(ScR, coeff, Length, Width, ARtable, mu,
                                                                           straindrop)
        MMO = (1 / c) * (np.log10(straindrop * mu * Length ** 2 * Width) - d)
        M = np.array([MMO, MAR, MRLD, MRA]);        M = np.round(M * 100) / 100
        dM1 = np.array([dMMO, ar_coeff[2], dMRLD, dMRA]);        dM1 = np.round(dM1 * 100) / 100
        M = M[~np.isnan(M)];         dM1 = dM1[~np.isnan(dM1)];         Mmean = np.mean(M)
        M1 = mag - Mmean
        # Conflation of Maximum magnitudes
        sdmagg=sdmag
        if mag < Mmean:
            flag_mobs = 1
            if abs(mag - Mmean) < Zeta:
                sdmag = sdmag
            elif abs(mag - Mmean) > Zeta:
                sdmag = np.mean(dM1) + Khi * abs(mag - Mmean)
        else:
            flag_mobs = 1
            if abs(mag - Mmean) < Zeta:
                sdmag = sdmag
            elif abs(mag - Mmean) > Zeta:
                print('Warning: Please consider revising the geometry parameters.')

        M = [MMO, MAR, MRLD, MRA, mag];        dM = [dMMO, ar_coeff[2], dMRLD, dMRA, sdmag]
        dM = np.round(np.array(dM) * 100) / 100  # Round to two decimal places
        # Correctly assigning a row in a NumPy array
        new_row = np.concatenate(([Length / 1000, LAR / 1000], M, dM)).reshape(1, -1)
        # Initialize M_dM_Lengths_forLOGoutput or stack the new row on top if it already exists
        if 'M_dM_Lengths_forLOGoutput' not in locals() or M_dM_Lengths_forLOGoutput.size == 0:
            M_dM_Lengths_forLOGoutput = new_row
        else:
            M_dM_Lengths_forLOGoutput = np.vstack((M_dM_Lengths_forLOGoutput, new_row))
        # Remove NaN values from M and dM
        M = [value for value in M if not np.isnan(value)]
        dM = [value for value in dM if not np.isnan(value)]
        # Convert M and dM to NumPy arrays if they are not already
        M = np.array(M)
        dM = np.array(dM)
        # Determine the minimum and maximum values considering M and dM
        min_val = np.floor(np.min(M - dM));
        max_val = np.ceil(np.max(M + dM))
        # Create x_range_of_mag with a step of 0.01
        step = 0.01
        x_range_of_mag = np.arange(min_val, max_val + step, step)
        ########################################### Including the n_sigma to truncate magnitudes:
        # Calculate a probability density function by a normal distribution for each M
        pdf_magnitudes = []
        for k in range(len(M)):
            pdf_magnitudes.append(norm.pdf(x_range_of_mag, M[k], dM[k]))
        pdf_magnitudes = np.array(pdf_magnitudes)
        pdf_magnitudes /= np.tile(np.max(pdf_magnitudes, axis=1)[:, np.newaxis], (1, pdf_magnitudes.shape[1]))
        # calcuate the summed distribution if LAR is greater or equal than Length then LAR is not used
        if LAR >= Length:
            pdf_magnitudes = np.delete(pdf_magnitudes, 1, axis=0)
        summed_pdf_magnitudes = np.sum(pdf_magnitudes, axis=0)
        ############## Calculate conflated distribution using the previously defined conflate_pdfs function
        conflated = conflate_pdfs(x_range_of_mag, pdf_magnitudes)



        if flag_mobs==1:
            pp= norm.pdf(x_range_of_mag, mag, sdmagg)

            pdf_magnitudes[-1] = pdf_magnitudes[-1] * (np.trapz(pp)/np.trapz(pdf_magnitudes[-1]))

        weighted_mean = np.average(x_range_of_mag, weights=summed_pdf_magnitudes)
        weighted_std = np.sqrt(np.average((x_range_of_mag - weighted_mean) ** 2, weights=summed_pdf_magnitudes))

        Mmax = weighted_mean
        sigma_Mmax = weighted_std
        pdf_test = norm.pdf(x_range_of_mag, Mmax, sigma_Mmax)
        # Round to the first decimal
        Mmax = round(Mmax * 10) / 10
        sigma_Mmax = round(sigma_Mmax * 10) / 10
        # Output the rounded values
        print(f"Mmax: {Mmax}, sigma_Mmax: {sigma_Mmax}")
        if 'output_Mmax_sigmaMmax' not in locals():
            output_Mmax_sigmaMmax = np.empty((0, 2), float)
        # Calculate the PDF using the normal distribution
        gauss_fit = norm.pdf(x_range_of_mag, loc=Mmax, scale=sigma_Mmax)
        ########################## PLOTS
        ii = fault_name
        namefig = 'Conflation_of_PDFs' + str(ii)
        # Set up the figure with your desired style
        # plt.style.use('dark_background')
        plt.style.use('default')
        plt.figure(figsize=(8, 6))
        count_pdf = 0
        plt.plot(x_range_of_mag, pdf_magnitudes[count_pdf, :], linewidth=1.2, linestyle='-', color='blue',
                 label='Label1')
        count_pdf += 1
        if LAR < Length:
            plt.plot(x_range_of_mag, pdf_magnitudes[count_pdf, :], linewidth=1.4, linestyle='-', color='lime',
                     label='Label2')
            count_pdf += 1
        plt.plot(x_range_of_mag, pdf_magnitudes[count_pdf, :], linewidth=1.0, linestyle='-', color='red',
                 label='Label3')
        count_pdf += 1
        plt.plot(x_range_of_mag, pdf_magnitudes[count_pdf, :], linewidth=1.2, linestyle='-', color='cyan',
                 label='Label4')
        count_pdf += 1
        if not np.isnan(mag):
            plt.plot(x_range_of_mag, pdf_magnitudes[count_pdf, :], linewidth=1.2, linestyle='-', color='magenta',
                     label='Label5')
        plt.plot(x_range_of_mag, summed_pdf_magnitudes, linewidth=1.2, linestyle='--', color='gray', label='Label6')
        plt.plot(x_range_of_mag, conflated, linewidth=2.0, linestyle='-', color='gold', label='Label7')
        aa = np.max(conflated)
        Mmax1 = x_range_of_mag[np.argmax(conflated)]

        # Make Mmax visible by changing the color to white
        plt.stem(Mmax1, np.max(conflated), linefmt='k-', markerfmt='ko', basefmt=' ')
        stem_proxy = mpatches.Patch(color='black', label='Mmax')
        # Check LAR and Length conditions and define legend entries
        if LAR < Length:
            if not np.isnan(mag):
                legendEntries = ['MMo'] + ['MAR'] + legends_Mw + ['MObs', 'SEM', 'CoP', 'Mmax']
            else:
                legendEntries = ['MMo'] + ['MAR'] + legends_Mw + ['SEM', 'CoP', 'Mmax']
        else:
            if not np.isnan(mag):
                legendEntries = ['MMo'] + legends_Mw + ['MObs', 'SEM', 'CoP', 'Mmax']
            else:
                legendEntries = ['MMo'] + legends_Mw + ['SEM', 'CoP'] + [] + ['Mmax']
        # Display the legend
        plt.legend(legendEntries)
        # plt.plot([Mmax1 - sigma_Mmax, Mmax1 + sigma_Mmax], [np.max(conflated), np.max(conflated)], linewidth=1.5,
        #          linestyle='-.', color='black')
        plt.fill_between(x_range_of_mag, 0, conflated, color='gold', alpha=0.25)
        plt.xlabel('Magnitude', fontsize=14, fontname='Times')
        plt.ylabel('Probability density function', fontsize=14, fontname='Times')
        plt.title(fault_name)
        # Plot settings:
        plt.ylim(0,)
        plt.xlim(5,)
        plt.subplots_adjust(left=0.075, right=0.98, top=0.96, bottom=0.1)  # Adjust margins of the plot
        # Save the figure

        if len(ProjFol) > 0:
            figure_name = "./" + ProjFol + "/Figures/" + namefig + '.pdf'
            if not os.path.exists("./" + ProjFol + "/Figures/"):
                # If it doesn't exist, create it
                os.makedirs("./" + ProjFol + "/Figures/")
        else:
            figure_name = "./" + 'output_files' + "/Figures/" + namefig + '.pdf'
            if not os.path.exists("./" + 'output_files' + "/Figures/"):
                # If it doesn't exist, create it
                os.makedirs("./" + 'output_files' + "/Figures/")

        # Save the figure in PDF format with higher DPI
        plt.savefig(figure_name, dpi=1200, bbox_inches='tight')  # Adjust the DPI value as needed
        # plt.title(ii)
        plt.show()  # Show the plot
        Mmax1 = round(Mmax1 * 10) / 10
        Mmax = Mmax1
        # Save the figure
        plt.savefig(namefig + '.pdf', dpi=1200, bbox_inches='tight')  # Adjust the DPI value as needed


        Mmax1 = round(Mmax1 * 10) / 10
        Mmax = Mmax1
        # Save the figure
        plt.savefig(namefig + '.png', dpi=600)  # Save the figure as a PNG file
        if LAR >= Length:
            L_forTmean = Length
        else:
            L_forTmean = LAR


        # Calculate Tmean
        Mmax=Mmax1
        # sigma_Mmax=0.5
        # Mmax=7.1
        output_Mmax_sigmaMmax = np.vstack([output_Mmax_sigmaMmax, [Mmax, sigma_Mmax]])
        Tmean = np.round(10 ** (d + c * Mmax) / (mu * V * L_forTmean * Width))  # Average recurrence time as defined in Field, 1999
        Tmean = np.round(Tmean * (1 / SCC))
        # Initialize Tep and varTep
        Tep = np.empty_like(Mmax)
        varTep = np.empty_like(Mmax)
        # Calculate Tep
        Tep = (10 ** (d + c * Mmax) / (mu * V * L_forTmean * Width) +
               (10 ** (d + c * Mmax) * c * np.log(10) / (mu * V * L_forTmean * Width)) * sigma_Mmax -
               10 ** (d + c * Mmax) / (mu * (V ** 2) * L_forTmean * Width) * dV)
        # Scale Tep by 1/SCC
        Tep *= (1 / SCC)
        # Calculate varTep
        varTep = ((10 ** (d + c * Mmax) * c * np.log(10) / (mu * V * L_forTmean * Width)) ** 2 * (sigma_Mmax ** 2) +
                  (10 ** (d + c * Mmax) / (mu * (V ** 2) * L_forTmean * Width)) ** 2 * (dV ** 2))
        # Scale varTep by (1/SCC)^2
        varTep *= (1 / SCC) ** 2
        # Calculate e
        e = np.sqrt(varTep)
        # Calculate alfa
        alfa = e / Tmean
        # Calculate MomentRate
        MomentRate = 10 ** (d + c * Mmax) / Tmean
        fault_data = {}
        # Assuming 'i' is an index variable
         # You can replace this with your actual index
        # Create a dictionary to store fault data
        fault_data = {
            "id": kk,  # Your 'i' index
            "Mmax": output_Mmax_sigmaMmax[kk - 1][0],  # Adjust the index, assuming 1-based indexing
            "sdMmax": output_Mmax_sigmaMmax[kk - 1][1],  # Adjust the index
            "Tmean": Tmean,
            "CV": alfa,
            "Telap": Telap,
            "MomentRate": MomentRate.tolist()  # Convert MomentRate to a Python list
        }
        kk=kk+1
        # Adding the outputs of the moment budget to the faults dictionary
        faults[fault_name].update(fault_data)

    return faults

def sactivityrate(faults, Fault_behaviour, w, bin, ProjFol):
    d = 9.1;  c = 1.5
    field_names = list(faults.keys())
    nfault = len(field_names)
    # Assuming 'faults' is your dictionary
    field_names = list(faults.keys())
    # Initializing arrays
    id = np.empty(len(field_names))
    mag = np.empty(len(field_names))
    sdmag = np.empty(len(field_names))
    Tmean = np.empty(len(field_names))
    alpha_val = np.empty(len(field_names))
    Telapsed = []
    Morate_input = np.empty(len(field_names))
    # ProbabilityOfOccurrence and Max_Telap aren't used in your MATLAB code
    fault_name = []
    idgr = np.empty(len(field_names))
    mt = np.empty(len(field_names))
    b = np.empty(len(field_names))
    # Iterate over each fault
    for i, field_name in enumerate(field_names):
        sub_struct = faults[field_name]
        id[i] = sub_struct.get('id', np.nan)
        mag[i] = sub_struct.get('Mmax', np.nan)
        sdmag[i] = sub_struct.get('sdMmax', np.nan)
        Morate_input[i] = sub_struct.get('MomentRate', np.nan)
        Tmean[i] = sub_struct.get('Tmean', np.nan)
        Telapsed.append(faults[field_name]['Telap'])
        alpha_val[i] = sub_struct.get('CV', np.nan)
        idgr[i] = sub_struct.get('id', np.nan)
        mt[i] = sub_struct.get('Mmin', np.nan)
        b[i] = sub_struct.get('b-value', np.nan)
        fault_name.append(field_name)

    Morate_fromTmean = 10 ** (c * mag + d) / Tmean
    Tmean_fromMorate = np.round(10 ** (c * mag + d) / Morate_input).astype(int)

    Morate = Morate_fromTmean

    for i in range(nfault):
        if Morate_fromTmean[i] != Morate_input[i]:
            print(
                f"Warning: Mo rate computed using M and Tmean for the fault # {i} is {Morate_fromTmean[i]:.4e}, different from Mo rate given in the input {Morate_input[i]:.4e}")
    Hbpt = np.zeros(nfault)  # Initialize the Hbpt array
    Hpois = np.zeros(nfault)
    for i in range(nfault):
        # sdmag[i]=0.5
        magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
        M = 10 ** (c * magnitude_range + d)
        pdf_mag = norm.pdf(magnitude_range, mag[i], sdmag[i])
        total_moment = np.sum(pdf_mag * M)
        ratio = Morate_input[i] / total_moment
        balanced_pdf_moment = ratio * pdf_mag
        CumRateMmin = np.sum(balanced_pdf_moment)
        Tm = 1 / CumRateMmin

        Telap = Telapsed[i]
        if not np.isnan(Telap):
            if Telap > 10 * Tm:
                Telap = 10 * Tm
                print(
                    f'Warning: Telap for fault id {id[i]} is forced to be equal to 10*Tm to avoid computational problems')
            alpha = alpha_val[i]

            scale = Tm / (alpha ** 2)
            # In Python: using scipy.stats.invgauss
            Hbpt_a1 = invgauss.cdf((Telap + w) / scale, mu=Tm / scale)
            Hbpt_a2 = invgauss.cdf((Telap) / scale, mu=Tm / scale)
            Hbpt[i] = (Hbpt_a1 - Hbpt_a2) / (1 - Hbpt_a2)
            if Hbpt[i] > 1:
                Hbpt[i] = 1

        # Assuming i, w, and Tm are defined earlier in your code
        Hpois[i] = 1 - np.exp(-1 * w * (1 / Tm))
        Hpois[Hpois > 1] = 1

    if Fault_behaviour == "Characteristic Gaussian" and Telapsed[i]:
        # bin=0.2
        CHGaussBPT(faults, c, d, ProjFol, fault_name, mag, sdmag, Tmean, Morate, id, nfault, w, Hbpt, bin)
    elif Fault_behaviour == "Characteristic Gaussian" and ~Telapsed[i]:
        CHGaussPoiss(faults, c, d, ProjFol, fault_name, mag, sdmag, Morate, id, nfault, w, Hpois, bin)
    elif Fault_behaviour == "Truncated Gutenberg Richter":
        TruncatedGR(faults, c, d, ProjFol, fault_name, mag, mt, Morate, id, nfault, bin, b)
    else:
        print("wrong case")
















