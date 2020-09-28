import numpy as np
from scipy import stats

# if running on Mac, sometimes the following is needed:
import matplotlib as mpl
mpl.use('TkAgg')

import matplotlib.pyplot as plt
import seaborn as sns

# keys are what the data is labeled as, values are what you'd like to see plots and such labeled as
pub_names = {'Trehalose':'Trehalose', 'Sugars':'Sugar Mix', 'Filterpaper':'Filter Paper',
             'Dry': 'No Additives', 'GenTegra':'GenTegra','PCR_DNAStable':'DNAStable + PCR',
             'DNAStable':'DNAStable', 'Beads_DNAStable':'Mag-Bind + DNAStable',
            'ETH_Trehalose':'ETH Trehalose', 'ETH_DNA_pure': 'ETH No Additives', 'ETH_Magnetic_NP':'ETH Magnetic NP',
             'ETH_DNA Stable':'ETH DNAStable',
            'Imagene':'Imagene'}

def open_csv(name, time):
    """Given the name of the .csv file and a list of the timepoints (in seconds) associated with it,
    returns the array of data associated with each temperature (tp1data, tp2data...)
    """

    # load .csv and put placeholders for missing values
    data = np.genfromtxt(name, delimiter=',', missing_values=np.nan)

    # ln of relative concentration
    nrows, ncols = data.shape
    norm_factor = data[:, 0].copy()
    for col in range(ncols):
        data[:, col] = data[:, col] / norm_factor
    data = np.log(data)

    # Appending together
    number_time_points = len(time)
    if number_time_points != ncols:
        raise Exception('Number of time points is incorrect for this csv')

    temp65_data, temp75_data, temp85_data = [], [], []
    for col in range(ncols):
        time_val = time[col]

        if 'ETH' in name:
            # 65C
            row = 0
            data_val = data[row, col]
            if not np.isnan(data_val):
                temp65_data.append([time_val, data_val])
            # 75C
            row = 1
            data_val = data[row, col]
            if not np.isnan(data_val):
                temp75_data.append([time_val, data_val])
            # 85C
            row = 2
            data_val = data[row, col]
            if not np.isnan(data_val):
                temp85_data.append([time_val, data_val])

        else:
            # Rows 1-3
            for row in range(0, 3):
                data_val = data[row, col]
                if not np.isnan(data_val):
                    temp65_data.append([time_val, data_val])

            # Rows 3-6
            for row in range(3, 6):
                data_val = data[row, col]
                if not np.isnan(data_val):
                    temp75_data.append([time_val, data_val])

            # Rows 6-9
            for row in range(6, 9):
                data_val = data[row, col]
                if not np.isnan(data_val):
                    temp85_data.append([time_val, data_val])

    temp65_data = np.array(temp65_data)
    temp75_data = np.array(temp75_data)
    temp85_data = np.array(temp85_data)

    return temp65_data, temp75_data, temp85_data


def linreg_conc_time(data, name, plot_command):
    # slope = -k
    slope, intercept, r_value, p_value, std_err = stats.linregress(data)

    if plot_command == 'yes_plot':
        plt.plot(data[:, 0], data[:, 1], 'bo', label=name)
        plt.plot(data[:, 0], slope*data[:, 0] + intercept, 'b', label='fit')
        plt.legend()
        plt.show()

    # returns k value
    return abs(slope)


def linreg_lnk_temp(k65, k75, k85, temp, plot_command):
    k_val = [k65, k75, k85]
    k_val = np.array(k_val)
    lnk_val = np.log(k_val)

    kelvin = np.array(temp) + 273.15
    inv_kelvin = 1.0/kelvin

    # where slope*8.314 is Ea, intercept is ln(k0)
    slope, intercept, r_value, p_value, std_err = stats.linregress(inv_kelvin, lnk_val)

    if plot_command == 'yes_plot':
        # plots ln(k) vs 1/T
        plt.plot(inv_kelvin, lnk_val, 'go', label='lnk v 1/T')
        plt.plot(inv_kelvin, slope*inv_kelvin + intercept, 'g-', label='fit')
        plt.legend()
        plt.show()

    # returns k0, -Ea/R, temp in kelvin
    return np.exp(intercept), slope, kelvin


def extrap_Arrhenius(k0, slope, temp_kel): # data [k0, Ea, temp_k]
    extrap_k = (k0 * np.exp(slope/temp_kel))
    halflife = (np.log(2) / extrap_k)

    # converts halflife from sec to years
    halflife_years = halflife/(60*60*24*365)

    return halflife_years

def get_colors():
    """"Returns a list of colors to use in plots that overlay data"""
    # http: // tools.medialab.sciences - po.fr / iwanthue /
    #colors = ["#cd9fba", "#86d64e","#8f44c8", "#d4c24d", "#7877c8", "#d5553a", "#8bd9a8", "#cb4a91", "#598140",
    #          "#4b2e5f", "#bb8f5e", "#75adbe", "#7e3737", "#414437"]
    colors = {'Trehalose': '#cd9fba', 'Sugars': '#86d64e', 'Filterpaper': '#8f44c8',
              'Dry': '#d4c24d', 'GenTegra': '#7877c8', 'PCR_DNAStable': '#d5553a',
              'DNAStable': '#8bd9a8', 'Beads_DNAStable': '#598140',
              'ETH_Trehalose': '#cb4a91', 'ETH_DNA_pure': '#bb8f5e', 'ETH_Magnetic_NP': '#4b2e5f',
              'ETH_DNA Stable': '#75adbe',
              'Imagene': '#7e3737'}
    return colors

def all_half_life_plots_together(list_of_treatment_files):
    """Given a list of csv file names (treatment_name.csv), returns
    a plot of all treatments and their half lives"""

    #time = [0, 191700, 340200, 513000, 1208700, 2417400]
    temp = [65, 75, 85]
    temp_new = np.arange(-15, 95, 10)
    temp_extra = temp_new + 273.15

    # makes plot pretty
    sns.set(style="white")
    sns.set_context("talk")
    fig, ax = plt.subplots(figsize=(5,8.5)) #(3.5, 5) is a nice skinny size #(12,5) for fat
    # http://tools.medialab.sciences-po.fr/iwanthue/

    colors = get_colors()

    # gets and plots data
    for i, f in enumerate(list_of_data_files):

        # adjusts to the fact that ETH and Imagene had different time points
        if 'ETH' in f:
            time = [0, 604800, 1209600, 1814400, 2419200, 3024000]
        elif 'Imagene' in f:
            time = [0, 191700, 513000, 1208700, 2417400, 3022200]
        else:
            time = [0, 191700, 340200, 513000, 1208700, 2417400]

        treatment = f.split('.')[0]  # gets the string that comes before the '.' in the filename
        temp65_dataset, temp75_dataset, temp85_dataset = open_csv(f, time)

        # gets the abs(slope) of data, AKA abs(k)
        # no longer shows plot unless specified 'yes_plot'
        reg65 = linreg_conc_time(temp65_dataset, '65', 'no_plot')
        reg75 = linreg_conc_time(temp75_dataset, '75', 'no_plot')
        reg85 = linreg_conc_time(temp85_dataset, '85', 'no_plot')

        # no longer shows plot unless specified 'yes_plot'
        k0, slope, temp_k = linreg_lnk_temp(reg65, reg75, reg85, temp, 'no_plot')

        # to account for DNA strand length
        if 'ETH' in f:
            #scale = 150.0
            scale = 1.0
        else:
            #scale = 310.0  # Illumina Truseq Nano kit with two 8bp adapters + 25N region
            scale = 310.0/150.0

        k0 *= scale

        k_extra = extrap_Arrhenius(k0, slope, temp_extra)  # k0 = A, used to call global var A

        half_lives = extrap_Arrhenius(k0, slope, temp_k)

        # plots our three data points (temp vs half_life)
        treatments_to_ignore = ['Imagene', 'ETH_DNAStable', 'ETH_DNA_pure', 'Dry']
        # treatments_to_analyze = ['ETH_Trehalose', 'Trehalose', 'ETH_Magnetic_NP', 'Sugars']

        if treatment not in treatments_to_ignore:
            if "ETH" not in treatment:
                ax.scatter(temp, half_lives, color=colors[treatment])

        # other than Imagene, plots extrapolated half life for other temps
        #if treatment != 'Imagene' and treatment != 'ETH_DNAStable' and treatment != 'ETH_DNA_pure':
        if treatment not in treatments_to_ignore and "ETH" not in treatment:
            if 'ETH' in treatment:
                ax.plot(temp_new, k_extra, color=colors[treatment], label=pub_names[treatment], linestyle='-.')
            if 'ETH' not in treatment:
                ax.plot(temp_new, k_extra, color=colors[treatment], label=pub_names[treatment], linestyle='-')


    # exceptions (data from other papers)
    ax.scatter(13, 800*(158.0/150), marker="^", c='#1E90FF', edgecolor= 'k', s=170, zorder=10, label="Moa Bone")  # DNA in moa bone (Allentoft et al. 2012)
    ax.scatter(25, 100*(158.0/150), marker="s", c="#1E90FF", edgecolor= 'k', s=160, zorder=10, label="Calcium Chloride")  # DNA over CaCl2 (Bonnet et al. 2010)
    ax.scatter(20, 80*(158.0/150), marker="X", c="#1E90FF", edgecolor= 'k', s=170, zorder=10, label="Silica")  # DNA in silica (Grass et al. 2010)
    ax.scatter(20, 0.8*(158.0/150), marker="p", c="#1E90FF", edgecolor= 'k', s=210, zorder=10, label="In Solution")  # DNA in solution (Lindahl & Nyberg 1972)

    # sets up rest of plot
    ax.set_yscale('log')
    ax.set_ylabel('Half-Life (Years)', fontsize='x-large')
    ax.set_xlabel('Temperature (C)', fontsize='x-large')
    ax.set_title('Half-Life', fontsize='x-large')
    sns.despine()

    # for legend in plot
    lgd = ax.legend(loc='upper right', frameon=True, facecolor="inherit", edgecolor="inherit",
              shadow=True , fontsize='medium', bbox_to_anchor=(2.6,1.02))

    # plt.show()  # NOTE- you cannot save an image when showing, you can either do one or the other
    plt.tight_layout()

    fig.savefig('skinny_UW_Half_Lives_per_nt.png', bbox_extra_artists=(lgd,),
                bbox_inches='tight', dpi=1000)

    plt.close()

def convert_sec_to_yrs(k):
    return k/(60*60*24*365)

<<<<<<< HEAD
=======

>>>>>>> cc4430cdb87119b7ee16bb739c916148980573ca
# roughly ordered via half-life projection
list_of_data_files = ['Dry.csv', 'Trehalose.csv', 'DNAStable.csv', 'Beads_DNAStable.csv',
                      'PCR_DNAStable.csv', 'Sugars.csv','Filterpaper.csv','GenTegra.csv',
                      'ETH_DNA_pure.csv', 'ETH_Trehalose.csv', 'ETH_DNAStable.csv', 'ETH_Magnetic_NP.csv',
                      'Imagene.csv']


all_half_life_plots_together(list_of_data_files)
