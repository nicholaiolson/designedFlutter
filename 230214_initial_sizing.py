# -*- coding: utf-8 -*-
"""
Nicholai Olson
230214

Purpose: Quick Initial Sizing Study
"""

import numpy as np
import dfFuncs as df
import matplotlib.pyplot as plt

#%% AEROSCOUT PARAMS
#https://www.horizonhobby.com/on/demandware.static/-/Sites-horizon-master/default/dw142e10ba/Manuals/HBZ380001_MANUAL_EN.pdf
c_aeroScout = 8*0.0254 #[m]
b_aeroScout = 42*0.0254 #[m]
S_aeroScout = c_aeroScout * b_aeroScout #[m^2]
W_aeroScout = 3.4*4.45 #[]
AR_aeroScout = b_aeroScout/c_aeroScout

W_S_aeroScout = W_aeroScout/S_aeroScout
# W_S_grob102 = 450/12.4*9.81

### TALON GT REBEL PARAMS
# https://www.getfpv.com/zohd-talon-gt-rebel-1000mm-wingspan-v-tail-bepp-fpv-aircraft-pnp.html?utm_source=google&utm_medium=cpc&utm_campaign=DM+-+B+-+PMax+-+Shop+-+SM+-+ALL&utm_content=pmax_x&utm_keyword=&utm_matchtype=&campaign_id=17881616054&network=x&device=c&gclid=CjwKCAiAmJGgBhAZEiwA1JZolllTh-SuyNVSnISFhzLcDT7XU86ao3K556Q4df9RbQoxfApfu0RnIBoCiOcQAvD_BwE

b_talonGT = 1 #[m]
S_talonGT = 0.14 #[m^2] c_talonGT * b_talonGT #[m^2]
c_talonGT =S_talonGT/b_talonGT
W_S_talonGT = 9.28*9.81#[N/m^2]
W_talonGT= W_S_talonGT*S_talonGT #[N]
AR_talonGT = b_talonGT/c_talonGT

### MUGIN-3 HTAIL PARAMS
# https://www.muginuav.com/product/mugin-3-3600m-h-tail-uav-platform/
b_mugin = 3.6 # [m]
W_mugin = 28*4.45 #[N]
S_mugin = 1.33 #[m^2]
c_mugin = S_mugin/b_mugin
W_S_mugin = W_mugin/S_mugin
AR_mugin = b_mugin/c_mugin


#%% BUILD SENSITIVITY ARRAYS
# wing loading / stall speed (takeoff performance), total wing span (I want it to fit into the Fit), and low frequencies.

npoints = 50

#Wing Loading 
minus = 0.8
plus =2.3
W_S = np.linspace(W_S_aeroScout*minus,W_S_aeroScout*plus,npoints)

#Weight
W_nominal = 9
minus = .67
plus = 1.33
M = 5.5
W = M*9.807

#REference Area
S = W/W_S

#Stall Speed
CLmax = 1.45
VS = (W/(CLmax*S))**(1/2)


#Chord
C_nominal = 7*.0254
minus = 0.67
plus = 1.2
C = .125#5*.0254
# C = S/b

#Span 
b = S/C

#Aspect Ratio
AR_nominal = 30
minus = 0.67
plus = 1.5
AR = b/C #np.linspace(AR_nominal*minus, nominal*plus, npoints)


plt.plot(W_S, b, label = 'bref')
plt.plot(W_S, S, label = 'Sref')
# plt.plot(W_S, C/.0254, label = 'Cref')
# plt.plot(W_S, AR, label = 'AR')
plt.plot([W_S_aeroScout*0.5, W_S_aeroScout*1.5], [AR_aeroScout, AR_aeroScout], label = "AeroScout AR")
plt.plot([W_S_aeroScout*0.5, W_S_aeroScout*1.5], [S_aeroScout, S_aeroScout], label = "AeroScout Sref")
plt.plot([W_S_aeroScout, W_S_aeroScout], [0, 30], label = "AeroScout W/S")

plt.grid()
plt.legend()

plt.figure(11)
plt.plot(W_S, VS, label = "Stalling Speed")
plt.plot(W_S, b, label = "Span")




#%% IMPORT AIRFOILS
clarky = {}
clarky_67P = {}
clarky['points_o'] = df.readAirfoil(r'C:\Users\nicho\designedFlutter\clarky.dat',1)
clarky_67P['points_o'] = df.readAirfoil(r'C:\Users\nicho\designedFlutter\clarky_67p.dat',1)




#Build up dictionary of desired wing sections to test - NEEDS TO KNOW CHORD
sections = []
sections.append({'description': 'Single Ply, 3K, 2 x 2 Twill Weave, ClarkY', 
                 'thickness':0.012*.0254, 
                 'airfoil': clarky['points_o']})
sections.append({'description': 'Single Ply, 3K, 2 x 2 Twill Weave, ClarkY 67%Thick', 
                 'thickness':0.012*.0254, 
                 'airfoil': clarky_67P['points_o']})

sections.append({'description': 'Dual Ply, 3K, 2 x 2 Twill Weave, ClarkY', 
                 'thickness':2*0.012*.0254, 
                 'airfoil': clarky['points_o']})
sections.append({'description': 'Dual Ply, 3K, 2 x 2 Twill Weave, ClarkY 67%Thick', 
                 'thickness':2*0.012*.0254, 
                 'airfoil': clarky_67P['points_o']})


for section in sections: 
    section['geometric properties'] = df.constantThickApproxAirfoilProps(section['airfoil'], C, section['thickness'])


#%% MODE ESTIMATES
#Mode shopaes for constant section cantilever beam, per Page 77 of BAH
E = 2.316634e+11 #Per https://s3.amazonaws.com/cdn.fibreglast.com/downloads/00101.pdf 
sigma_ult = 4205.8e6 #Per https://s3.amazonaws.com/cdn.fibreglast.com/downloads/00101.pdf 
rho = 1670 #kg/m^3
massScaling = 2
nz = 9


for section in sections:
    m = massScaling*b/2*section['geometric properties']['A']*rho #Mass of semispan
    omega1 = (0.597)**2*(3.1415**2/(b/2)**2)*(E*section['geometric properties']['Ixx_c']/m)**(1/2)/6.28
    omega2 = (1.49)**2*(3.1415**2/(b/2)**2)*(E*section['geometric properties']['Ixx_c']/m)**(1/2)/6.28
    
    omega1_halfEI_m = (0.597)**2*(3.1415**2/(b/2)**2)*(E*section['geometric properties']['Ixx_c']/(m*2))**(1/2)/6.28
    
    omega1_pointMass = 1/(2*3.1415)*(3*E*section['geometric properties']['Ixx_c']/(m*(b/2)**3))**(1/2)
    
    #Loads and stress
    Mx = b/2*.4*nz*W/2
    sigma = Mx*((np.max(section['airfoil'][:,1])-np.min(section['airfoil'][:,1]))/2)/section['geometric properties']['Ixx_c']
    
    plt.figure(7)
    plt.plot(b,sigma, label = section['description'])
    
    plt.figure(8)
    plt.plot(b,omega1, label = section['description'])

    plt.figure(9)
    plt.plot(W_S,omega1, label = section['description'])
    

plt.figure(7)
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.plot(b,sigma_ult*np.ones(b.shape),label='Ultimate Tensile Strength')
plt.legend()
plt.grid()
plt.title(str(round(W/9.807,1)) +' kg Airplane, ' + str(round(C,3)) + ' m chord')
plt.xlabel('Span, b [m]')
plt.ylabel('Bending Stress [Pa]')
plt.legend(loc='upper left')


plt.figure(8)
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.legend()
plt.grid()
plt.title(str(round(W/9.807,1)) +' kg Airplane, ' + str(round(C,3)) + ' m chord')
plt.xlabel('Span, b [m]')
plt.ylabel('First Wing Bending Mode [Hz]')
plt.legend(loc='upper left')
plt.figure(9)
plt.legend()
plt.grid()


plt.figure(11)
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.plot(W_S,omega1, label = 'First Wing Bending Mode, ' + section['description'])
plt.plot(W_S,omega1_halfEI_m, label = 'First Wing Bending Mode, ' + section['description'] + ' Half EI/m')
plt.plot(W_S,omega1_pointMass, label = 'First Wing Bending Mode, ' + section['description'] + ' Point Mass')
plt.plot([W_S_aeroScout, W_S_aeroScout], [0, max(omega1)], label = "AeroScout W/S")
plt.plot([W_S_talonGT, W_S_talonGT], [0, max(omega1)], label = "Talon GT W/S")
plt.plot([W_S_mugin, W_S_mugin], [0, max(omega1)], label = "Mugin-3 W/S")
# plt.plot([W_S_grob102, W_S_grob102], [0, 15], label = "AeroScout W/S")

plt.title(str(round(W/9.807,1)) +' kg Airplane, ' + str(round(C,3)) + ' m chord')
plt.xlabel('Wing Loading, W/S [N/m^2]')
plt.ylabel('Equivalent Stalling Speed [m/s], Wing Span [m], Frequency [Hz]')
plt.grid()
plt.legend(loc='upper left')

#%% LOADS & STRESS
#Compute wing root bending assuming CP at 40% span



#%% WING WEIGHT ESTIMATE