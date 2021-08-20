# Get the bands from a quantum espresso pw.x "bands" calculation. Must be used with verbosity = 'high'
# The output has a shape nk x (nband + 1), where nk is the number of k points and nband is the number of bands.
# The first column is the distance along the k path and the rest of the columns are the energies. 
# The resulting output can then easily be plotted with matplotlib. For example, if the output file is "bands.dat"
# then the bands can be plotted with:

# bands = np.loadtxt('bands.dat')
# ks = bands[:,0]
# for i in range(1,bands.shape[1]):
#     plt.plot(ks, bands[:,i], 'k')

#the x axis of the plot is the distance along the k path, which is stored in the first column of the output

import numpy as np
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "number of Kohn-Sham states" in line:
        nband = int(line.split(' ')[-1])

    if "number of k points" in line:
        nk = int(line.split('=')[1].split()[0])
        
    if "End of band structure calculation" in line:
        start_index = i+1
        break

bands = np.zeros((nk, nband+1))

k_index = -1
k_last = np.zeros(3)
band_index = 1
for line in lines[start_index:]:
    if "Writing output data file" in line:
        break
    if "k =" in line:
        k_index += 1
        band_index = 1
        kline = line.split('=')[1].split()
        kx = float(kline[0])
        ky = float(kline[1])
        kz = float(kline[2])
        kdist = np.sqrt((kx - k_last[0])**2 + (ky - k_last[1])**2 + (kz - k_last[2])**2)
        if k_index == 0:
            bands[k_index, 0] = kdist
        else:
            bands[k_index, 0] = kdist + bands[k_index - 1, 0]
        k_last[0] = kx
        k_last[1] = ky
        k_last[2] = kz
    else:
        for energy in line.split():
            bands[k_index, band_index] = float(energy)
            band_index += 1
        
np.savetxt(outfile, bands)
