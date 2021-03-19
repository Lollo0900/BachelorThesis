import pylab
import numpy
import random
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sympy import symbols, solve

cos_theta_cm,fhi_cm,cos_thetax_cm,fhix_cm,pa_zero_cm,pa_x_cm,pa_y_cm,pa_z_cm,pa_zero_lab,pa_x_lab,pa_y_lab,pa_z_lab,px_zero_cm,px_x_cm,px_y_cm,px_z_cm,px_zero_lab,px_x_lab,px_y_lab,px_z_lab,modpx,fl,flt=pylab.loadtxt('Generator\output.txt',unpack=True)
#plot istogrammi

plt.figure(1)
plt.title(r'$\tau$ $cos(\theta)$ in C.M.')
plt.xlabel(r'cos($\theta$)')
plt.ylabel('Occurrencies')
plt.hist(cos_theta_cm, bins = 20)
plt.savefig("Figure_1.pdf")


plt.figure(2)
plt.title(r'$\phi$ $\tau$ in C.M.')
plt.xlabel(r'$\phi$ [rad]')
plt.ylabel('Occurrencies')
plt.hist(fhi_cm, bins = 10)
plt.savefig("Figure_2.pdf")

plt.figure(3)
plt.title(r"$p_{\tau}$-x component in C.M.")
plt.xlabel(r'$p_{\tau}^x$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_x_cm, bins = 100)
plt.savefig("Figure_3.pdf")

plt.figure(4)
plt.title(r"$p_{\tau}$-x component in LAB")
plt.xlabel(r'$p_{\tau}^x$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_x_lab, bins = 100)
plt.savefig("Figure_4.pdf")

plt.figure(5)
plt.title(r"$p_{\tau}$-y component in C.M.")
plt.xlabel(r'$p_{\tau}^y$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_y_cm, bins = 100)
plt.savefig("Figure_5.pdf")

plt.figure(6)
plt.title(r"$p_{\tau}$-y component in LAB")
plt.xlabel(r'$p_{\tau}^y$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_y_lab, bins = 100)
plt.savefig("Figure_6.pdf")

plt.figure(7)
plt.title(r"$p_{\tau}$-z component in C.M.")
plt.xlabel(r'$p_{\tau}^z$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_z_cm, bins = 100)
plt.savefig("Figure_7.pdf")

plt.figure(8)
plt.title(r"$p_{\tau}$-z component in LAB")
plt.xlabel(r'$p_{\tau}^z$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(pa_z_lab, bins = 100)
plt.savefig("Figure_8.pdf")

plt.figure(9)
plt.title(r'$\Lambda^0$ $cos(\theta)$ in C.M.')
plt.ylabel('Occurrencies')
plt.hist(cos_thetax_cm, bins = 10)
plt.savefig("Figure_9.pdf")

plt.figure(10)
plt.title(r"$\phi$ $\Lambda^0$ in $\tau$ frame")
plt.xlabel(r'$\phi$ [rad]')
plt.ylabel('Occurrencies')
plt.hist(fhix_cm, bins = 10)
plt.savefig("Figure_10.pdf")

plt.figure(11)
plt.title(r"$p_{\Lambda^0}$-x component in $\tau$ frame")
plt.xlabel(r'$p_{\Lambda^0}^x$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_x_cm, bins = 100)
plt.savefig("Figure_11.pdf")

plt.figure(12)
plt.title(r"$p_{\Lambda^0}-x$ component in LAB frame")
plt.xlabel(r'$p_{\Lambda^0}^x$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_x_lab, bins = 100)
plt.savefig("Figure_12.pdf")

plt.figure(13)
plt.title(r"$p_{\Lambda^0}-y$ component in $\tau$ frame")
plt.xlabel(r'$p_{\Lambda^0}^y$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_y_cm, bins = 100)
plt.savefig("Figure_13.pdf")

plt.figure(14)
plt.title(r"$p_{\Lambda^0}-y$ component in LAB frame")
plt.xlabel(r'$p_{\Lambda^0}^y$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_y_lab, bins = 100)
plt.savefig("Figure_14.pdf")

plt.figure(15)
plt.title(r"$p_{\Lambda^0}-z$ component in $\tau$ frame")
plt.xlabel(r'$p_{\Lambda^0}^z$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_z_cm, bins = 100)
plt.savefig("Figure_15.pdf")

plt.figure(16)
plt.title(r"$p_{\Lambda^0}-z$ component in LAB frame")
plt.xlabel(r'$p_{\Lambda^0}^z$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(px_z_lab, bins = 100)
plt.savefig("Figure_16.pdf")

plt.figure(17)
plt.title(r'$\Lambda^0$ 3-momenta module in C.M. frame')
plt.xlabel(r'$p_{\Lambda^0}$ [GeV/c]')
plt.ylabel('Occurrencies')
plt.hist(modpx, bins = 50, range=(1.5,max(modpx)))
plt.savefig("Figure_17.pdf")
plt.xticks(numpy.arange(1.5, max(modpx), 0.15))

plt.figure(18)
plt.title(r'Flight length $L_{\Lambda^0}$ in LAB')
plt.xlabel(r'$L_{\Lambda^0}$ [cm]')
plt.ylabel('Occurrencies')
plt.hist(fl, bins = 500)
plt.savefig("Figure_18.pdf")

plt.figure(19)
plt.title(r'Flight length $L_{\tau}$ in LAB')
plt.xlabel(r'$L_{\tau}$ [cm]')
plt.ylabel('Occurrencies')
plt.hist(flt, bins = 500)
plt.savefig("Figure_19.pdf")

plt.show()
