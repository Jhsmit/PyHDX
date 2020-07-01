from math import exp, log10

# Here are the parameters measured in Bai et al. (1993)
# The parameters for D and E are based on the work of Mori and al. (1997)
# and measured by and calculated in the functions acid and base
para={
    "A": [0.00, 0.00, 0.00, 0.00],
    "C": [-0.54, -0.46, 0.62, 0.55],
    "F": [-0.52, -0.43, -0.235859464059171, 0.0631315866300978],
    "G": [-0.22, 0.218176047120386, -0.03, 0.17],
    "I": [-0.91, -0.59, -0.73, -0.23],
    "K": [-0.56, -0.29, -0.04, 0.12],
    "L": [-0.57, -0.13, -0.576252727721677, -0.21],
    "M": [-0.64, -0.28, -0.00895484265292644, 0.11],
    "N": [-0.58, -0.13, 0.49, 0.38],
    "P": [9999, -0.194773472023435, 9999, -0.24],
    "Q": [-0.47, -0.27, 0.06, 0.20],
    "R": [-0.59, -0.32, 0.0767122542818456, 0.22],
    "S": [-0.437992277698594, -0.388518934646472, 0.37, 0.299550285605933],
    "T": [-0.79, -0.448073125742265, -0.0662579798400606, 0.20],
    "V": [-0.739022273362575, -0.30, -0.701934483299758, -0.14],
    "W": [-0.40, -0.44, -0.41, -0.11],
    "Y": [-0.41, -0.37, -0.27, 0.05],
}

rho_Nterm_acid = -1.32
rho_Nterm_base = 1.62

lamb_Cterm_acid=0.96
lamb_Cterm_base=-1.80

pKD = 15.05
R = 1.987

ka = 10**(1.62)/60
kb = 10**(10.18)/60
kw = 10**(-1.5)/60

Ea = 14000
Eb = 17000
Ew = 19000


lamb_cterm_acid=0.96

lamb_cterm_base=-1.80




def get_D(pH):
    return 10**(-pH)


def get_OD(pH):
    return 10**(pH-pKD)


def get_l_DTR(temperature):
    return (1./temperature-1/293.)/R


def get_pK_his(temperature):
    Ea_his = 7500
    return -log10(10**(-7.42)*exp(-Ea_his*(1/(R*temperature)-1/(278*R))))


def get_pK_asp(temperature):
    Ea_asp = 1000
    return -log10(10**(-4.48)*exp(-Ea_asp*(1/(R*temperature)-1/(278*R))))


def get_pK_glu(temperature):
    Ea_glu = 1083
    return -log10(10**(-4.93)*exp(-Ea_glu*(1/(R*temperature)-1/(278*R))))


def get_Fta(temperature):
    return exp(-Ea*get_l_DTR(temperature))


def get_Ftb(temperature):
    return exp(-Eb*get_l_DTR(temperature))


def get_Ftw(temperature):
    return exp(-Ew*get_l_DTR(temperature))
