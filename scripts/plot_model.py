
import config
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats



ds = numpy.logspace(-5, 0, base=10, num=100)


a = 0.1
b=0.5
dnds = (1-a) + a*(1-numpy.exp(-b*ds))/(ds*b)

fig = plt.figure(figsize = (4, 4))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax = plt.subplot2grid((1, 1), (0, 0))

ax.plot(ds, dnds)

ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)


fig.subplots_adjust(hspace=0.25, wspace=0.25)
fig_name = "%stest_model.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()