import tables as tb
import numpy as np


class Particle():
    def __init__(self, name, idnumber, ADCcount, TDCcount, grid_i, grid_j, pressure, energy):
        self.name = name
        self.idnumber = idnumber
        self.ADCcount = ADCcount
        self.TDCcount = TDCcount
        self.grid_i = grid_i
        self.grid_j = grid_j
        self.pressure = pressure
        self.energy = energy


h5file = tb.open_file("tutorial1.h5", mode="w", title="Test file")
group = h5file.create_group("/", 'detector', 'Detector information')
table = h5file.create_table(group, 'readout', Particle, "Readout example")

particle = table.row
for i in range(10):
    particle['name'] = 'Particle: %6d' % (i)
    particle['TDCcount'] = i % 256
    particle['ADCcount'] = (i * 256) % (1 << 16)
    particle['grid_i'] = i
    particle['grid_j'] = 10 - i
    particle['pressure'] = float(i*i)
    particle['energy'] = float(particle['pressure'] ** 4)
    particle['idnumber'] = i * (2 ** 34)
    # Insert a new particle record
    particle.append()

table.flush()
