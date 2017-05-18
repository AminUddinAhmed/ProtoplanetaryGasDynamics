# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 12:02:17 2017

@author: The Users
"""
import time
import numpy as np
import matplotlib.pyplot as plt

start_time = time.time()

#declaring constants
universalGasConstant = 8.31441
gravitationalConstant = 6.67259 * 10**-11
astronomicalUnit = 1.496 * (10**11)
earthMass = 5.976 * 10**24
sunMass = 1.989 * (10**30)
coreRadius = 0
coreDensity = 5500.0 / (2**9)
surfaceDensity = 0
semiMajorAxis = 5.2 * astronomicalUnit
coreMass = earthMass * 10**-5
gasMeanMolecularMass = 2.3 * 10**-3
sunSurfaceTemperature = 5800
sunRadius = 6.96 * 10**8
equilibriumTemperature = sunSurfaceTemperature * \
    (sunRadius / (2 * semiMajorAxis))**(1 / 2)
print('equilibrium temperature at', semiMajorAxis,
      'metres is', equilibriumTemperature, 'kelvin')
muOverRT = gasMeanMolecularMass/(equilibriumTemperature * universalGasConstant)
print('muOverRT is', muOverRT)

#creating arrays
initialEnvelopeDensity = np.array([])
finalEnvelopeDensity = np.array([])
massSolutions = np.array([])
densitySolutions = np.array([])

earths = 34
points = 280

for mass in range(0, earths):    
    coreMass = coreMass * 1.5
    
    for log in range(0, points):
                
        initialDensityLog = log - 125
        initialDensity = 1.09**initialDensityLog 
        print(mass,'/', earths, log,'/',points, 'initial density is', initialDensity)
        #surfaceDensity = initialDensity * 10**initialDensityLog
        initialEnvelopeDensity = np.append(initialEnvelopeDensity, initialDensity)
    
        
        #print('coremass is', coreMass)
        coreRadius = (coreMass / ((4 / 3) * np.pi * coreDensity))**(1 / 3)
        #print('coreradius is', coreRadius)
        
        hillSphere = (semiMajorAxis) * (((coreMass) / (3 * sunMass))**(1 / 3))
        #print('hillsphere is', hillSphere)#THIS IS THE HILLSPHERE
        
        numberOfSlices = 2000
        deltaRadius = hillSphere / numberOfSlices
        #print('slice length is', deltaRadius)
        
        
        totalMass = coreMass
        
        integrationRange = 2000
        densityAtRadius = initialDensity
        deltaDensity = 0
        deltaMass = 0
        pressureAtRadius = 0
        distanceAxis = np.array([])
        densityAxis = np.array([])
        massAxis = np.array([])
        hillSphereAxis = np.array([])
        distanceFromCOG = coreRadius
        pressureAxis = np.array([])
        massAxis = np.array([])
        #print ('totalmass initially is', totalMass)
        
        
        for slice in range(0,integrationRange):   
          
                        
            hillSphere = (semiMajorAxis) * (((totalMass) / (3 * sunMass))**(1 / 3))
            hillSphereAxis = np.append(hillSphereAxis, hillSphere)
            deltaRadius = hillSphere / numberOfSlices
            
            distanceFromCOG = coreRadius+(slice*deltaRadius)
            distanceAxis = np.append(distanceAxis, distanceFromCOG)    
            
            #Runge-Kutta algorithm to integrate
        
            k1r = distanceFromCOG
            k1rho = densityAtRadius
            k1m  = totalMass
            k1density = -(((gravitationalConstant*k1m*k1rho*muOverRT))/(k1r**2))
            k1mass = 4*np.pi*k1rho*(((k1r)**2))
        
            k2r = distanceFromCOG + 0.5*deltaRadius
            k2rho = densityAtRadius + 0.5*k1density*deltaRadius
            k2m  = totalMass + 0.5*k1mass*deltaRadius
            k2density = -(((gravitationalConstant*k2m*k2rho*muOverRT))/(k2r**2))
            k2mass = 4*np.pi*k2rho*(((k2r)**2))
            
            k3r = distanceFromCOG + 0.5*deltaRadius
            k3rho = densityAtRadius + 0.5*k2density*deltaRadius
            k3m  = totalMass + 0.5*k2mass*deltaRadius
            k3density = -(((gravitationalConstant*k3m*k3rho*muOverRT))/(k3r**2))
            k3mass = 4*np.pi*k3rho*(((k3r)**2))
            
            k4r = distanceFromCOG + deltaRadius
            k4rho = densityAtRadius + k3density*deltaRadius
            k4m  = totalMass + k3mass*deltaRadius
            k4density = -(((gravitationalConstant*k4m*k4rho*muOverRT))/(k4r**2))
            k4mass = 4*np.pi*k4rho*(((k4r)**2))    
            
            
            nextDensity = densityAtRadius + (1/6.0)*(k1density + 2*k2density + 2*k3density + k4density)*deltaRadius
            densityAtRadius = nextDensity
            densityAxis = np.append(densityAxis, nextDensity)
            nextMass = totalMass + (1/6.0)*(k1mass + 2*k2mass + 2*k3mass + k4mass)*deltaRadius
            totalMass = nextMass
            massAxis = np.append(massAxis, nextMass)        
            
            
            
        
        finalEnvelopeDensity = np.append(finalEnvelopeDensity, densityAxis[-1])
        print(densityAxis[-1])
        if 8*10**-7 < densityAxis[-1] < 1.2*10**-6:
            densitySolutions = np.append(densitySolutions, initialDensity)
            massSolutions = np.append(massSolutions, coreMass)
        '''
        print("new total mass is", totalMass)
        print("envelope mass is", totalMass - coreMass)
        print (densityAxis)
        
        print(densityAxis.shape)
        print(distanceAxis.shape)
        #print (massAxis)
        print(hillSphereAxis)
        
        #print(distanceAxis)
        
        pressureAxis = muOverRT * densityAxis
        
        plt.plot(distanceAxis, densityAxis)
        plt.xlabel('Radius')
        plt.ylabel('Density')
        plt.title('Density as function of radius')
        plt.grid(True)
        plt.show()
        
        
        plt.plot(distanceAxis, hillSphereAxis)
        plt.xlabel('Radius')
        plt.ylabel('hill radius')
        plt.title('hillsphere as function of radius')
        plt.grid(True)
        plt.show()
        
        plt.plot(distanceAxis, massAxis)
        plt.xlabel('Radius')
        plt.ylabel('Mass')
        plt.title('Mass as function of radius')
        plt.grid(True)
        plt.show()
        '''
    
'''
plt.loglog(distanceAxis, pressureAxis)
plt.xlabel('Radius')
plt.ylabel('Pressure')
plt.title('Pressure as function of radius')
plt.grid(True)
plt.show()

'''

#testing array sizes
print(initialEnvelopeDensity.shape)
print(initialEnvelopeDensity)


print(finalEnvelopeDensity.shape)
print(finalEnvelopeDensity)


#graphs below

plt.semilogx(initialEnvelopeDensity, finalEnvelopeDensity)
plt.xlabel('surface density')
plt.ylabel('outer density')
plt.ylim(0, 0.00005)
#plt.grid(True)
plt.show()

plt.loglog(initialEnvelopeDensity, finalEnvelopeDensity)
plt.xlabel('surface density')
plt.ylabel('outer density')
plt.grid(True)
plt.show()

print(densitySolutions)
print(massSolutions)


plt.scatter(massSolutions, densitySolutions)
axes = plt.gca()
axes.set_xscale('log')
axes.set_yscale('log')
#plt.xlim(10**19, 10**26)
plt.ylim(10**-6, 10**6)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
#plt.grid(True)
plt.show()



plt.scatter(massSolutions, densitySolutions)
axes = plt.gca()
axes.set_xscale('log')
axes.set_yscale('log')
plt.ylim(1, 1000000)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()
'''
plt.scatter(massSolutions, densitySolutions)
axes = plt.gca()
axes.set_xscale('log')
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()

plt.scatter(massSolutions, densitySolutions)
axes = plt.gca()
axes.set_yscale('log')
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()

plt.loglog(massSolutions, densitySolutions)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()

plt.semilogx(massSolutions, densitySolutions)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()

plt.semilogy(massSolutions, densitySolutions)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()

plt.scatter(massSolutions, densitySolutions)
plt.xlabel('mass')
plt.ylabel('density')
plt.title('critical')
plt.grid(True)
plt.show()
'''


print("--- %s seconds ---" % (time.time() - start_time))
