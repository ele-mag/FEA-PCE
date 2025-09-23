import sys

import numpy as np
import pandas as pd
import skrf as rf
from pylab import *
from skrf import Network, Frequency, NetworkSet
def read_sdB_2_re_parameters(file_path):

    data = pd.read_csv(file_path, comment='!', delimiter='\s+',
                       names=['Frequency', 'dbS11', 'angS11', 'dbS21', 'angS21', 'dbS12', 'angS12', 'dbS22', 'angS22'])

    data = data.apply(pd.to_numeric, errors='coerce')

    data = data.dropna()

    def db_to_mag(db):
        return 10 ** (db / 20)

    for col in ['S11', 'S21', 'S12', 'S22']:
        data[f'ang{col}'] = np.radians(data[f'ang{col}'])
        data[f'db{col}'] = db_to_mag(data[f'db{col}'])

    for col in ['S11', 'S21', 'S12', 'S22']:
        data[f'complex{col}'] = data[f'db{col}'] * np.exp(1j * data[f'ang{col}'])

    frequency = rf.Frequency.from_f(data['Frequency'].values, unit='GHz')

    s_parameters_complex = np.array([
        [data['complexS11'].values, data['complexS21'].values],
        [data['complexS12'].values, data['complexS22'].values]]).transpose((2, 0, 1))

    network = rf.Network(frequency=frequency, s=s_parameters_complex, z0=50)

    return network
    # for col in ['S11', 'S21', 'S12', 'S22']:
    #     data[f'real{col}'] = data[f'complex{col}'].apply(lambda x: x.real)
    #     data[f'imag{col}'] = data[f'complex{col}'].apply(lambda x: x.imag)
    # return data[['Frequency', 'realS11', 'imagS11', 'realS21', 'imagS21', 'realS12', 'imagS12', 'realS22', 'imagS22']]
jianhe_netw = read_sdB_2_re_parameters("Micr.s2p")
jianhe_netw.plot_s_db()

plt.show()
