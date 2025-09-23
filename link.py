import sys
import pandas as pd
import skrf as rf
from pylab import *
import define_lianlu
from skrf import Network, Frequency, NetworkSet
import bonding_scan2snp
from write_s import fuction


f_min = 8
f_max = 12
freq = Frequency(f_min, f_max, 500, 'ghz')
def load_and_resample_network(path, name, frequency):
    net = rf.Network(path)
    net.name = name
    net.resample(frequency)
    return net

network_files = {
    'Phase': 'TR_SYS/s_python/phase/HMC543ALC4B_sparm_25C_1_0.s2p',
    'Sw_A': 'TR_SYS/s_python/SW/25C/ADRF5025_RF1_SEL_CTRL_3.3_25C.s3p',
    'Filter_TA': r'TR_SYS/s_python/filter/filter.s2p',
    'Pa_TA': 'TR_SYS/s_python/PA1/HMC1082LP4E_S-par_includes_EVB_111213.s2p',
    'ATT_T': 'TR_SYS/s_python/ATT/HMC424A-DIE_direct probed sparameters_0.0dB.s2p',
    'Pa_TB': 'TR_SYS/s_python/PA2/HMC952ALP5GE_SParameters.s2p',
    'Filter_TB': 'TR_SYS/s_python/filter/filter.s2p',
    'Sw_B': 'TR_SYS/s_python/SW/25C/ADRF5025_RF1_SEL_CTRL_3.3_25C.s3p',
    'dianqiao': 'TR_SYS/s_python/3db/3db.s4p',
    'Filter_TC': r'TR_SYS/s_python/filter/filter.s2p',
    'Filter_R': 'TR_SYS/s_python/filter/filter.s2p',
    'ATT_RA': 'TR_SYS/s_python/ATT/HMC424A-DIE_direct probed sparameters_0.0dB.s2p',
    'LNA1': 'TR_SYS/s_python/LNA/HMC963LC4 De-embedded.s2p',
    'LNA2': 'TR_SYS/s_python/LNA/HMC963LC4 De-embedded.s2p',
    'ATT_RB': 'TR_SYS/s_python/ATT/HMC424A-DIE_direct probed sparameters_0.0dB.s2p',
    'LNA_B': 'TR_SYS/s_python/LNA/HMC963LC4 De-embedded.s2p',
    'ssma': 'TR_SYS/s_python/ssma.s2p',
}

for name, path in network_files.items():
    net = load_and_resample_network(path, name, freq)
    globals()[name] = net

network_files_line = {
    "l1": 'bonding1.s2p',
    "l2": 'bonding2.s2p',
    "l3": 'bonding3.s2p',
    "l4": 'bonding4.s2p',
    "l5": 'bonding5.s2p',
    "l6": 'bonding6.s2p',
    "l7": 'bonding7.s2p',
    "l8": 'bonding8.s2p',
    "l9": 'bonding9.s2p',
    "l10": 'bonding10.s2p',
    "l11": 'bonding11.s2p',
    "l12": 'bonding12.s2p',
}

def TR_R():

    port1 = rf.Circuit.Port(freq, name="port1", z0=50)
    port2 = rf.Circuit.Port(freq, name="port2", z0=50)
    port3 = rf.Circuit.Port(freq, name="port3", z0=50)
    port1.resample(freq)
    port2.resample(freq)
    port3.resample(freq)
    ground1 = rf.Circuit.Ground(freq, name='ground1')
    jiedi_dianzu = rf.media.DefinedGammaZ0(frequency=freq, z0=50)
    R1 = jiedi_dianzu.resistor(50, name='R1')

    for index, (name, path) in enumerate(network_files_line.items()):
        a = a_values[index]
        net = read_sdB_2_re_parameters(path, name=name, beishu=a)

        globals()[name] = net

    connexions = [

        [(port1, 0), (ssma1, 0)],
        [(ssma1, 1), (l1, 0)],
        [(l1, 1), (Sw_A, 0)],
        [(Sw_A, 1), (l2, 0)],
        [(l2, 1), (Filter_TA, 0)],
        [(Filter_TA, 1), (l3, 0)],
        [(l3, 1), (Pa_TA, 0)],
        [(Pa_TA, 1), (l4, 0)],
        [(l4, 1), (Filter_TB, 0)],
        [(Filter_TB, 1), (GF, 0)],
        [(GF, 1), (l6, 0)],
        [(GF, 1), (l10, 0)],

        [(l6, 1), (Filter_TC, 0)],
        [(Filter_TC, 1), (L7, 0)],
        [(L7, 1), (LNA1, 0)],
        [(LNA1, 1), (l8, 0)],
        [(l8, 1), (Sw_B, 1)],
        [(l10, 1), (LNA2, 0)],
        [(LNA2, 1), (l11, 0)],
        [(l11, 1), (lim, 0)],
        [(lim, 1), (l12, 0)],
        [(l12, 1), (Sw_B, 2)],
        [(Sw_B, 0), (l9, 0)],
        [(l9, 0), (ssma2, 1)],
        [(ssma2, 1), (port2, 1)],

    ]

    cir = rf.Circuit(connexions)
    net = cir.network

    value = net.s_db[:, 0, 1]
    # value = net.s_db[:, 1, 0]
    print(freq.f, value)
    return freq.f, value

