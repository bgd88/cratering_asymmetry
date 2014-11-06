import pymc as pm
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# Define key parameters, Values from Zahnle 1998
vp = 13.070  # Orbital velocity of Jupiter, [km/s]
R_j =  71398.0 # jovian radii [km]
G = 6.67*10**(-11) # Universal Gravitational Constant 

moon = 'Ganymede'

v_orb = {'Io' : 17.3, 
         'Europa' : 13.7,
         'Ganymede' : 10.9,
         'Callisto' : 8.2}  # Orbital velocity of Jupiter's Moons, [km/s]

R_sat = {'Io' : 1820,
         'Europe' : 1570,
         'Ganymede' : 2630,
         'Callisto' : 2400} # Radius [km]

a_sat = {'Io' : 5.9*R_j,
         'Europa' : 9.4*R_j,
         'Ganymede' : 15.0*R_j,
         'Callisto' : 26.4*R_j} # Semimajor axis of orbit [km]

g = {'Io' : 1.81,
     'Europa' : 1.31,
     'Ganymede' : 1.43,
     'Callisto' : 1.25} # Surface Gravity [m/s2]

def v_esc_eval(g,r):
    return np.sqrt(2*g*r*1000)/1000 # Calculates surface velocity 

x1 = pm.Uniform('x1', lower=0.0, upper=1.0) # Uniformly distributed random number
x2 = pm.Uniform('x2', lower=0.0, upper=1.0) # Uniformly distributed random number
x3 = pm.Uniform('x3', lower=0.0, upper=1.0) # Uniformly distributed random number
x4 = pm.Uniform('x4', lower=0.0, upper=1.0) # Uniformly distributed random number
x5 = pm.Uniform('x5', lower=0.0, upper=1.0) # Uniformly distributed random number
x6 = pm.Uniform('x6', lower=0.0, upper=1.0) # Uniformly distributed random number

def Vinf_eval(vp=vp,x1=x1): 
    return 0.3*vp*(1-0.45* np.log(1/(x1+0.2)-0.832))

def Uinf_eval(vp=vp,x1=x1,v_orb=v_orb[moon]):
    return Vinf_eval(vp,x1)/v_orb

#inclinations
def inclination(x2=x2):
    return np.arccos(1-2*x2)

#eccentricitys
#WARNING, THERE IS AN ERROR IN ZAHNLE 2001
def eccentricity(Uinf=0,x3=x3):
    UU_inf = Uinf**2
    return np.sqrt(1+x3*UU_inf*(UU_inf+2)) 

def azimuth_eval(x4=x4):
    return np.pi*(1-2*x4)

def theta_eval(x5=x5):
    return np.arccos(np.sqrt(x5))

#In the absence of gravitational focusing theta is equal to psi
# We have defined az as the chi in Zahnale
def colat_eval(az,psi,omega):
    return np.arccos(np.cos(psi)*np.cos(omega)+np.sin(psi)*np.sin(omega)*np.cos(az))

def lambda_eval(Uy,Ux):
    return(-Uy/Ux)

def omega_eval(Uz,U):
    return np.arccos(Uz/U)

def delta_eval(az,psi,omega,colat):
    sinD = np.sin(psi)*np.sin(az)/np.sin(colat)
    cosD = (np.cos(psi)-np.cos(colat)*np.cos(omega))/(np.sin(colat)*np.sin(omega))
    return np.arctan(sinD/cosD)                                              
                           
def v_enc_eval(U,v_orb):
    return U*v_orb

################################
#Gravitational Focusing
def psi_eval(v_enc,v_orb,x5=x5):
    nu = 2*(v_enc/v_orb)**2
    qp = (np.sqrt(1+x5*nu*(nu+2))-1)/nu
    ep = 1 + nu*qp
    return np.arccos(-1/ep)-np.arccos((qp*(1+ep)-1)/ep)

################################
#periapse
def periapse(Uinf,e):
    UU_inf = Uinf**2
    return (e-1)/UU_inf

def U_eval(e,q,i):
    return np.sqrt(3-(1-e)/q-2*np.cos(i)*np.sqrt(q*(1+e)))

def Ux_eval(e,q):
    return np.sqrt(2-(1-e)/q-q*(1+e))

def Uy_eval(e,q,i):
    return np.cos(i)*np.sqrt(q*(1+e))-1

def Uz_eval(e,q,i):
    return np.sin(i)*np.sqrt(q*(1+e))

def impact_probabilities(i,U,Ux,U_inf,R_sat,a_sat,v_orb,v_esc,R_j):
    Pstar = ( R_sat/a_sat )**2 *( 1+(v_esc/(v_orb*U))**2 )*(U/Ux)/(np.pi*np.sin(i))
    w = (2+ U_inf**2)/(2 + (R_j/a_sat)*U_inf**2 )
    return w*Pstar

#@pm.deterministic
def impact_event(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,vp=vp,v_orb=v_orb[moon]):
    v_esc = v_esc_eval(g=g[moon],r=R_sat[moon])
    Uinf = Uinf_eval(vp,x1,v_orb)
    e = eccentricity(Uinf,x3)
    q = periapse(Uinf,e)
    i = inclination(x2)
    Ux = Ux_eval(e,q)
    Uy = Uy_eval(e,q,i)
    Uz = Uz_eval(e,q,i)
    U = U_eval(e,q,i)
    v_enc = v_enc_eval(U,v_orb)
    omega = omega_eval(Uz,U) 
    lamb = lambda_eval(Uy,Ux)
    az = azimuth_eval(x4)
    psi = psi_eval(v_enc,v_esc,x5)
    #psi = theta_eval(x5) #psi i= theta for no gravitational forcing
    colat = colat_eval(az,psi,omega)
    delta =  delta_eval(az,psi,omega,colat)
    lat = np.pi/2-colat
    lon = delta + lamb
    beta = np.arccos(np.cos(lat)*np.sin(lon))
    P = impact_probabilities(i,U,Ux,Uinf,R_sat=R_sat[moon],a_sat=a_sat[moon],v_orb=v_orb,v_esc=v_esc,R_j=R_j)
    #loc = locals()
    return np.mod(beta,2*np.pi), P

num = 50000
b = np.empty((num,1))
p = np.empty_like(b)
for ii in np.arange(num):
    b[ii], p[ii] = impact_event(x1.random(),x2.random(),x3.random(),x4.random(),x5.random(),vp=vp,v_orb=v_orb[moon])

fig = plt.figure(figsize=(10,10))
h1, bins, patches = plt.hist(b*180/np.pi,50,normed=True,alpha=0.5,label='Apex Angle')
h2, bins, patches = plt.hist(b*180/np.pi,50,weights=p,normed=True,alpha=0.5,label='Weighted by Impact Prob')
plt.legend(loc='upper right')
plt.show()

area_weights = np.cos(bins[0:-1]*np.pi/180)-np.cos(bins[1:]*np.pi/180)

plt.bar(bins[0:-1],h2/area_weights)
plt.show()
