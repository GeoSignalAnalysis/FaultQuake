import math
import numpy
import numpy as np
import os
from scipy.stats import norm
from src.FaultQuake_functions import kin2coeff, coeff2mag, conflate_pdfs
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import invgauss

import statsmodels.api as sm


# def momentbudget():
# ShearModulusFromInputFile = ShearModulusFromInput
# StrainDropFromInputFile = StrainDropFromInput
# if ShearModulusFromInputFile is None:
#     ShearModulusFromInputFile = 0
# if StrainDropFromInputFile is None:
#     StrainDropFromInputFile = 0
# warning_off()
# if ShearModulus is None and ShearModulusFromInputFile != 1:
#     print("Warning: Shear modulus - You are using the default value 3*10^10 Pa")
#     ShearModulus = 3
# if StrainDrop is None and StrainDropFromInputFile != 1:
#     print("Warning: Strain Drop - You are using the default value 3*10^-5")
#     StrainDrop = 3
# if n_sigma is None or np.isnan(n_sigma):
#     print("Warning: You are using untruncated Gaussian distributions for magnitudes")
#     n_sigma_trunc_magnitudes = 0
# else:
#     n_sigma_trunc_magnitudes = n_sigma
# outputfile = os.path.join('./output_files', f'{outputname}_MB.txt')
# with open(outputfile, 'w') as fidout:
#     fidout.write('id Mmax sdMmax Tmean CV Telap Mo-rate name\n')
# temporary_pdf_storage = open('./output_files/temporary_pdf_storage.txt', 'w')
# print('RUNNING...')

# for item_name in faults:
#     print(f"Item Name: {item_name}")
# # mu_opt, straindrop_opt = ShearModulus, StrainDrop
# # if ShearModulusFromInputFile != 1:
#     mu_opt = ShearModulus * 1e10
# # if StrainDropFromInputFile != 1:
#     straindrop_opt = StrainDrop * 1e-5

# with open(imfile, 'r') as fid2:
#     raw = fid2.read()
# DATA = json.loads(raw)
# fault_name = list(DATA.keys())


# def momentbudget():
#
#     faults = self.browse_file()
#
#     for item_name in faults:
#         print(f"Item Name: {item_name}")
#         mu_opt = ShearModulus * 1e10
#         straindrop_opt = StrainDrop * 1e-5


def momentbudget(faults,ProjFol):
    # aaaa=1
    # for item_name in faults:
    #     print(f"Item Name: {item_name}")
    #     mu_opt = ShearModulus * 1e10
    #     straindrop_opt = StrainDrop * 1e-5
    #     # Use mu_opt and straindrop_opt as needed in your calculations.
    # for key in faults.items():
    #     # print(f"Variable Name: {variable_name}, Variable Value: {variable_value}")
    #    if faults[key]['ShearModulus'] = NAN
    #     math.isnan(value):
    #     if variable_name == 'ShearModulus':
    #         mu_opt = variable_value * 1e10
    #     elif variable_name == 'StrainDrop':
    #         straindrop_opt = variable_value * 1e-5

    # for key, value in faults.items():
    #     if value == 'ShearModulus' or value == 'StrainDrop':
    #         if math.isnan(faults[key]):
    #             ShearModulus = 3 * 1e10
    #         else:
    #             ShearModulus = faults['ShearModulus'] * 1e10
    #
    #         if math.isnan(faults[key]):
    #             StrainDrop = 3 * 1e-5
    #         else:
    #             StrainDrop = faults['StrainDrop'] * 3 * 1e-5
    #

    for fault, values in faults.items():
        if 'ShearModulus' in values and values['ShearModulus'] == 'NaN':
            print(f"Fault {fault} has a 'NaN' ShearModulus value.")
            values['ShearModulus'] = 3 * 1e10
        else:
            values['ShearModulus'] = values['ShearModulus'] * 1e10

        if 'StrainDrop' in values and values['StrainDrop'] == 'NaN':
            print(f"Fault {fault} has a 'NaN' ShearModulus value.")
            values['StrainDrop'] = 3 * 1e-5
        else:
            values['StrainDrop'] = values['StrainDrop'] * 1e-5
    kk=1
    for fault_name in faults:
        ScR = faults[fault_name]['ScR']
        # Length = np.append(Length, faults[fault_name]['Length'])
        Length = faults[fault_name]['Length']
        Dip = faults[fault_name]['Dip']
        Seismogenic_thickness = faults[fault_name]['Seismogenic_Thickness']
        Slipmin = faults[fault_name]['SRmin']
        Slipmax = faults[fault_name]['SRmax']
        mag = faults[fault_name]['Mobs']
        sdmag = faults[fault_name]['sdMobs']
        Last_eq_time = faults[fault_name]['Last_eq_time']
        SCC = faults[fault_name]['SCC']
        ShearModulusFromInputFile = faults[fault_name]['ShearModulus']
        StrainDropFromInputFile = faults[fault_name]['StrainDrop']
        yfc = faults[fault_name]['year_for_calculations']

        Telap = yfc - Last_eq_time  # Perform element-wise subtraction
        mu = ShearModulusFromInputFile
        straindrop = StrainDropFromInputFile

        Dip_radians = math.radians(Dip)  # Convert degrees to radians
        sine_Dip = math.sin(Dip_radians)

        Length = Length * 1000
        Width = (Seismogenic_thickness * 1000) / sine_Dip

        V = (Slipmin + Slipmax) / 2000
        dV = V - (Slipmin / 1000)

        # DEFINITION OF COEFFICIENTS
        # coefficients by Hanks & Kanamori, 1985
        d = 9.1
        c = 1.5
        dMMO = 0.3  # standard deviation of the magnitude calculated by the seismic moment definition

        # Calculation of the coefficient of the imported scale-relationship
        coeff, ARtable = kin2coeff(ScR)

        MRLD, MRA, dMRLD, dMRA, MAR, ar_coeff, LAR, legends_Mw = coeff2mag(ScR, coeff, Length, Width, ARtable, mu,
                                                                           straindrop)

        MMO = (1 / c) * (np.log10(straindrop * mu * Length ** 2 * Width) - d)

        M = np.array([MMO, MAR, MRLD, MRA])
        M = np.round(M * 100) / 100
        dM1 = np.array([dMMO, ar_coeff[2], dMRLD, dMRA])
        dM1 = np.round(dM1 * 100) / 100
        M = M[~np.isnan(M)]
        dM1 = dM1[~np.isnan(dM1)]
        Mmean = np.mean(M)
        M1 = mag - Mmean

        # Conflation of Maximum magnitudes

        if mag < Mmean:
            if abs(mag - Mmean) < 0.5:
                sdmag = sdmag
            elif abs(mag - Mmean) > 0.5:
                sdmag = np.mean(dM1) + 0.2 * abs(sdmag - np.mean(dM1))
        else:
            if abs(mag - Mmean) > 1:
                print('Warning: Please consider revising the geometry parameters.')

        M = [MMO, MAR, MRLD, MRA, mag]
        dM = [dMMO, ar_coeff[2], dMRLD, dMRA, sdmag]

        dM = np.round(np.array(dM) * 100) / 100  # Round to two decimal places

        # Combine M, dM, and Length/LAR for LOG output
        M_dM_Lengths_forLOGoutput = [Length / 1000, LAR / 1000, *M, *dM]

        # Remove NaN values from M and dM
        M = [value for value in M if not np.isnan(value)]
        dM = [value for value in dM if not np.isnan(value)]

        # M = np.array([MMO, MAR, MRLD, MRA, mag])
        # dM = np.array([dMMO, ar_coeff[2], dMRLD, dMRA, sdmag])
        # dM = np.round(dM * 100) / 100
        #
        # M = M[~np.isnan(M)]
        # dM = dM[~np.isnan(dM)]



        # Convert M and dM to NumPy arrays if they are not already
        M = np.array(M)
        dM = np.array(dM)

        # Determine the minimum and maximum values considering M and dM
        min_val = np.floor(np.min(M - dM))
        max_val = np.ceil(np.max(M + dM))

        # Create x_range_of_mag with a step of 0.01
        step = 0.01
        x_range_of_mag = np.arange(min_val, max_val + step, step)

        ########################################### Including the n_sigma to truncate magnitudes:

        # # Calculate Mmax
        # xmag = None
        # y = None
        # sumy = None
        # muhat = None
        # sigmahat = None
        # sdm = None
        # # x_range_of_mag = []
        #

        # if n_sigma_trunc_magnitudes == 0:
        #     x_range_of_mag = np.arange(np.floor(np.min(M - dM)), np.ceil(np.max(M + dM)) + 0.01, 0.01)
        # else:
        #     x_range_of_mag = np.arange(
        #         np.floor(np.min(M - n_sigma_trunc_magnitudes * dM)),
        #         np.ceil(np.max(M + n_sigma_trunc_magnitudes * dM)) + 0.01,
        #         0.01,
        #     )

        #################################################################################


        # Calculate a probability density function by a normal distribution for each M
        pdf_magnitudes = []

        for k in range(len(M)):
            pdf_magnitudes.append(norm.pdf(x_range_of_mag, M[k], dM[k]))


        ########################################### Including the n_sigma to truncate magnitudes:

        # pdf_magnitudes = calculate_pdf_magnitudes(M, dM, mag, n_sigma_trunc_magnitudes)
        # if n_sigma_trunc_magnitudes > 0:
        #     dist_trunc_mag = truncGaussDist(pdf_magnitudes, x_range_of_mag, M, dM, mag, n_sigma_trunc_magnitudes)
        #     pdf_magnitudes[0:dist_trunc_mag.shape[0], :] = dist_trunc_mag

        ###########################################################################################
        # Normalizing each value to of the pdfs to its maximum value, so each pdf weights as 1
        pdf_magnitudes = np.array(pdf_magnitudes)
        pdf_magnitudes /= np.tile(np.max(pdf_magnitudes, axis=1)[:, np.newaxis], (1, pdf_magnitudes.shape[1]))

        # calcuate the summed distribution if LAR is greater or equal than Length then LAR is not used


        summed_pdf_magnitudes = np.sum(pdf_magnitudes, axis=0)


        ############## Calculate conflated distribution using the previously defined conflate_pdfs function

        conflated = conflate_pdfs(x_range_of_mag, pdf_magnitudes)

        # Fit a normal distribution to the data in summed_pdf_magnitudes


        # Assuming x_range_of_mag is your data
        # summed_pdf_magnitudes is your weights or probabilities

        # Fit a normal distribution to the data
        Mmax = np.sum(x_range_of_mag * summed_pdf_magnitudes) / np.sum(summed_pdf_magnitudes)

        # Compute the weighted standard deviation
        weighted_diff = (x_range_of_mag - Mmax) ** 2 * summed_pdf_magnitudes
        sigma_Mmax = np.sqrt(np.sum(weighted_diff) / np.sum(summed_pdf_magnitudes))

        # Now, 'Mmax' contains the weighted mean, and 'sigma_Mmax' contains the weighted standard deviation
        norm_dist = norm(loc=Mmax, scale=sigma_Mmax)

        # Assuming x_range_of_mag is your range of values, calculate the PDF
        pdf_test = norm_dist.pdf(x_range_of_mag)


        # Round to the first decimal
        Mmax = round(Mmax * 10) / 10
        sigma_Mmax = round(sigma_Mmax * 10) / 10

        # Create a list to save the values
        # Initialize the output array if it's not already defined
        if 'output_Mmax_sigmaMmax' not in locals():
            output_Mmax_sigmaMmax = np.empty((0, 2), float)

        # Append [Mmax, sigma_Mmax] to the array as a new row
        output_Mmax_sigmaMmax = np.vstack([output_Mmax_sigmaMmax, [Mmax, sigma_Mmax]])

        # Calculate the PDF using the normal distribution
        gauss_fit = norm.pdf(x_range_of_mag, loc=Mmax, scale=sigma_Mmax)



        ########################## PLOTS


        ii = fault_name
        namefig = 'Conflation_of_PDFs' + str(ii)
        # plt.figure(ii)

        # Set up the figure with your desired style
        plt.style.use('dark_background')
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
        plt.stem(Mmax1, np.max(conflated), linefmt='w-', markerfmt='wo', basefmt=' ')
        stem_proxy = mpatches.Patch(color='white', label='Mmax')
        # Check LAR and Length conditions and define legend entries
        if LAR < Length:
            if not np.isnan(mag):
                legendEntries = ['MMo'] + ['MAR'] + legends_Mw + ['MObs', 'sEM', 'CoP', 'Mmax']
            else:
                legendEntries = ['MMo'] + ['MAR'] + legends_Mw + ['sEM', 'CoP', 'Mmax']
        else:
            if not np.isnan(mag):
                legendEntries = ['MMo'] + legends_Mw + ['MObs', 'sEM', 'CoP', 'Mmax']
            else:
                legendEntries = ['MMo'] + legends_Mw + ['sEM', 'CoP'] + [] + ['Mmax']
        # Display the legend
        plt.legend(legendEntries)
        # plt.plot([Mmax1 - sigma_Mmax, Mmax1 + sigma_Mmax], [np.max(conflated), np.max(conflated)], linewidth=1.5,
        #          linestyle='-.', color='black')
        plt.fill_between(x_range_of_mag, 0, conflated, color='gold', alpha=0.25)
        plt.xlabel('Magnitude', fontsize=14, fontname='Times New Roman')
        plt.ylabel('Probability density function', fontsize=14, fontname='Times New Roman')
        # Plot settings:
        plt.ylim(0,)
        # plt.xlim(5,)
        plt.subplots_adjust(left=0.075, right=0.98, top=0.96, bottom=0.1)  # Adjust margins of the plot
        # Save the figure
        if len(ProjFol)>0 :
            figure_name="./"+ProjFol+"/Figures/"+ namefig + '.eps'
            if not os.path.exists("./"+ProjFol+"/Figures/"):
                # If it doesn't exist, create it
                os.makedirs("./"+ProjFol+"/Figures/")
        else:
            figure_name="./"+'Outputs'+"/Figures/"+ namefig + '.eps'
            if not os.path.exists("./"+'Outputs'+"/Figures/"):
                # If it doesn't exist, create it
                os.makedirs("./"+'Outputs'+"/Figures/")



        plt.savefig(figure_name, dpi=800)  # Save the figure as a PNG file
        # plt.title(ii)
        plt.show()  # Show the plot
        # Save the figure
        plt.savefig(namefig + '.png', dpi=600)  # Save the figure as a PNG file
        if LAR >= Length:
            L_forTmean = Length
        else:
            L_forTmean = LAR
        # Calculate Tmean
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

    # Assuming all the required variables are defined and available
    # nfault, fixed_smr, Morate_fromTmean, Morate_input, Tmean_fromMorate, Tmean, etc.



    for i in range(nfault):
        if Morate_fromTmean[i] != Morate_input[i]:
            print(
                f"Warning: Mo rate computed using M and Tmean for the fault # {i} is {Morate_fromTmean[i]:.4e}, different from Mo rate given in the input {Morate_input[i]:.4e}")

    for i in range(nfault):
        magnitude_range = np.arange(mag[i] - sdmag[i], mag[i] + sdmag[i] + bin, bin)
        M = 10 ** (c * magnitude_range + d)
        pdf_mag = norm.pdf(magnitude_range, mag[i], sdmag[i])
        total_moment = np.sum(pdf_mag * M)
        ratio = Morate_input[i] / total_moment
        balanced_pdf_moment = ratio * pdf_mag
        CumRateMmin = np.sum(balanced_pdf_moment)
        Tm = 1 / CumRateMmin



        Hbpt = np.zeros(nfault)  # Initialize the Hbpt array
        Hpois = np.zeros(nfault)
        Telap = Telapsed[i]
        if not np.isnan(Telap):
            if Telap > 10 * Tm:
                Telap = 10 * Tm
                print(
                    f'Warning: Telap for fault id {id[i]} is forced to be equal to 10*Tm to avoid computational problems')
            alpha = alpha_val[i]
            scale = Tm / (alpha ** 2)

            # In MATLAB: cdf('inversegaussian', (Telap+w), Tm, (Tm/(alpha^2)))
            # In Python: using scipy.stats.invgauss
            Hbpt_a1 = invgauss.cdf((Telap + w) / scale, mu=Tm / scale)
            Hbpt_a2 = invgauss.cdf((Telap) / scale, mu=Tm / scale)
            Hbpt[i] = (Hbpt_a1 - Hbpt_a2) / (1 - Hbpt_a2)
            if Hbpt[i] > 1:
                Hbpt[i] = 1

        # Assuming i, w, and Tm are defined earlier in your code
        Hpois[i] = 1 - np.exp(-1 * w * (1 / Tm))
        Hpois[Hpois > 1] = 1



