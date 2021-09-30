import numpy
import matplotlib.pyplot as plt
import matplotlib.colors

pixels = numpy.load("pixels.npy")
plt.imshow(pixels, norm=matplotlib.colors.LogNorm(), cmap='hot')
plt.colorbar()
plt.show()
