import numpy as np
from numba import jit

@jit
def sigmaEff(v, threshold=0.683):
	v = np.sort(v)

	total = len(v)
	Max = int(threshold * total)

	start, stop, width = [], [], []

	i = 0
	while (i != len(v)-1):
		count = 0
		j = i
		while (j != len(v)-1 and count < Max):
			count +=1
			j += 1

		if (j != len(v)-1):
			start.append(v[i])
			stop.append(v[j])
			width.append(v[j] - v[i])

		i += 1

	npminwidth = np.array(width)

	minwidth = np.amin(npminwidth)
	pos = np.argmin(npminwidth)

	xmin = start[pos]
	xmax = stop[pos]

	return xmin, xmax, minwidth*0.5

# np.random.seed(2020)
# data = np.random.normal(0, 15, 1000)
# print(sigmaEff(data, 0.683))