import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 

carac = 'scccc'

H = np.linspace(0.1,0.99,100)
P_H = np.ones(100)

def after_vero(carac,H,P_H):
	carac_a = list(carac)
	N_c = 0
	N_s = 0
	for element in carac_a:
		if element == 's':
			N_s += 1
		else:
			N_c += 1

	vero = (H**(N_c))*((1-H)**N_s)

	after = vero*P_H

	P_obs = integrate.trapz(after,H)

	return after/float(P_obs)

after = after_vero(carac,H,P_H)
L = np.log(after)

first_prime = (L[1:]-L[:-1])/float(H[1]-H[0])
second_prime = (first_prime[1:]-first_prime[:-1])/float(H[1]-H[0])

ind = np.where(abs(first_prime) <= 0.15)

H_0 = float(H[ind])
sigma = float((-second_prime[ind])**(-0.5))

def normal(x,mu,sig):
	return np.exp(-((x-mu)**2)/(2*sigma**2))/((2*np.pi*sigma**2)**(0.5))

H_normal = np.linspace(0,1.5,100)

plt.figure()
plt.plot(H,after)
plt.plot(H_normal,normal(H_normal,H_0,sigma), '--', label = 'Gaussian')
plt.title(r'$H_0: {:.2f} \pm {:.2f}$'.format(H_0,sigma))
plt.xlabel(r'$H$')
plt.ylabel(r'$P(H|{obs})$')
plt.legend()
plt.savefig('coins.png')
