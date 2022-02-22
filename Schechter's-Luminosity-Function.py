import pandas as pd 
from matplotlib import pyplot as plt
import math as mt 
import numpy as np
from scipy.optimize import curve_fit

# Loading the data
df = pd.read_csv('SDSS.csv', encoding='ISO-8859-1')

# Magnitude color diagram before clearing data

# Defining some essential constants
H0 = 70        # km/s/Mpc (Hubble constant)
c = 3*10**5    # km/s (light velocity)

df['dist'] = df['redshift'] * c / H0 # Mpc
df['Mg'] = df['g'] - 5*np.log10(10**6*df['dist']) + 5
Mg = np.array(df['Mg'])
redshift = np.array(df['redshift'])
g = np.array(df['g'])
r = np.array(df['r'])
cor = np.array(g-r)
dist = np.array(df['dist'])

grafic_cond = np.where(((g-r)<1.0) & ((g-r)>0.0))


plt.scatter(Mg[grafic_cond], (g-r)[grafic_cond], s=1)
plt.xlabel('magnitude')
plt.ylabel('(g-r)')
plt.yticks(np.arange(0, 1.2, .2))
plt.title('color-magnitude diagram')
plt.show()

# Magnitude color diagram after clearing data
Mg_sample = Mg
cor_sample = cor

# selection criteria
restric = np.where( (Mg < -16) & (cor < 0.6) )

Mg = Mg[restric]
redshift = redshift[restric]
cor = cor[restric]

# Blue dots
grafic_cond_rest = np.where((cor<1.0) & (cor>0.0) )

# Red dots
grafic_cond = np.where(((g-r)<1.0) & ((g-r)>0.0))

# Gray dots

grafic_cond2 = np.where(((g-r)<1.0) & ((g-r)>0.0) & (Mg_sample>-16))

plt.scatter(Mg_sample[grafic_cond], cor_sample[grafic_cond], s=0.08, color='red')
plt.scatter(Mg_sample[grafic_cond2], cor_sample[grafic_cond2], s=0.2, color='gray')
plt.scatter(Mg[grafic_cond_rest], cor[grafic_cond_rest], s=0.5, color='blue')

plt.plot([-16, -16], [0.0,1.0], color='black')
plt.plot([-22, -10], [0.6, 0.6], color='black')
plt.xlabel('magnitude')
plt.ylabel('(g-r)')
plt.yticks(np.arange(0, 1.2, .2))
plt.gca().invert_xaxis()
plt.title('color-magnitude diagram')
plt.show()

# Preparing the code to calculate Schechter's Luminosity Function
bin = (np.max(Mg) - np.min(Mg))/15 # intervalo
i = np.min(Mg)
num_gal = []
Mg_mean = []
while i < np.max(Mg):
    i_plus = i+bin
    cond = np.where( (Mg>=i) & (Mg<i_plus) )
    num_gal.append(len(Mg[cond]))
    Mg_mean.append(np.mean(Mg[cond]))
    i=i+bin

    
num_gal = np.array(num_gal)
Mg_mean = np.array(Mg_mean)


dmax = np.max(dist) # dmax em Mpc
Omega = mt.pi # ângulo sólido
V = Omega/3 * dmax**3

phi = num_gal/(V*bin)

# Fit of the Schechter's Luminosity Function

#define function
def schechterM(magnitude, phiStar, alpha, MStar): 
    """Schechter luminosity function by magnitudes.""" 
    MStarMinM = 0.4 * (MStar - magnitude) 
    return np.log10(0.4 * np.log(10.0) * phiStar 
            * 10.0**(MStarMinM * (alpha + 1.0)) * np.exp(-10.0**MStarMinM)) 

#fit
params, cov = curve_fit(schechterM, Mg_mean, np.log10(phi), p0 = [0.1, 0, -20])

phiStar = params[0]
alpha   = params[1]
MStar   = params[2]

print(f'phiStar = {phiStar}')
print(f'alpha   = {alpha}')
print(f'MStar   = {MStar}')
LStar = 10**((-MStar+4.77)/2.5)
print(f'LStar   = {LStar}')
#evaluate func using fitted parameters
xfit = np.linspace(-21, -16,20)
yfit = schechterM(xfit, phiStar, alpha, MStar)

#plot data and fit
fig, ax = plt.subplots()

plt.plot(Mg_mean, np.log10(phi), 'o', label='data')
plt.plot(xfit, yfit, label='fit')

plt.gca().invert_xaxis()

plt.legend()

plt.xlabel('magnitude')
plt.ylabel('log($\phi$) [Mpc$^{-3}$ mag$^{-1}$]')
plt.title("Schechter's Luminosity Function")

plt.show()