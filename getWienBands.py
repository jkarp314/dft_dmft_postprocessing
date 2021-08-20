# Get the bands from a wien2k spaghetti calculation.
# The output has a shape nk x (nband + 1), where nk is the number of k points and nband is the number of bands.
# The first column is the distance along the k path and the rest of the columns are the energies. 
# The resulting output can then easily be plotted with matplotlib. For example, if the output file is "bands.dat"
# then the bands can be plotted with:

# bands = np.loadtxt('bands.dat')
# ks = bands[:,0]
# for i in range(1,bands.shape[1]):
#     plt.plot(ks, bands[:,i], 'k')

#the x axis of the plot is the distance along the k path, which is stored in the first column of the output

#example of running script:
# python3 getWienbands.py case.spaghetti_ene bands.dat
import numpy as np
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as f:
    lines = f.readlines()
nk = 0
for line in lines[1:]:
    if 'bandindex' in line:
        break
    nk += 1
nband = len(lines)//(nk+1)

bands = np.zeros((nk,nband+1))

k = 0
band_index = 1
for line in lines[1:]:
    if 'bandindex' in line:
        k = 0
        band_index += 1
        continue
    split_line = [x for x in line.split(' ') if x != ''] 
    bands[k, 0] = float(split_line[-2])
    bands[k, band_index] = float(split_line[-1])
    k += 1

np.savetxt(outfile, bands)
