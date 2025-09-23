import skrf as rf
def resample2snp(path, name, frequency):
    net = rf.Network(path)
    net.name = name
    net.resample(frequency)
    return net