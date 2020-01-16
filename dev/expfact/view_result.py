import matplotlib.pyplot as plt
import numpy as np


res, pf = np.genfromtxt('output/testoutput.pfact').T

plt.figure()
plt.scatter(res, pf)
plt.show()
