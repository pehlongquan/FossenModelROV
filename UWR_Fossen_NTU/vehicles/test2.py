import numpy as np
import math
import sys
import sys


def Smtrx(a):
    """
    S = Smtrx(a) computes the 3x3 vector skew-symmetric matrix S(a) = -S(a)'.
    The cross product satisfies: a x b = S(a)b. 
    """
 
    S = np.array([ 
        [ 0, -a[2], a[1] ],
        [ a[2],   0,     -a[0] ],
        [-a[1],   a[0],   0 ]  ])

    return S

def Hmtrx(r):
    """
    H = Hmtrx(r) computes the 6x6 system transformation matrix
    H = [eye(3)     S'
         zeros(3,3) eye(3) ]       Property: inv(H(r)) = H(-r)

    If r = r_bg is the vector from the CO to the CG, the model matrices in CO and
    CG are related by: M_CO = H(r_bg)' * M_CG * H(r_bg). Generalized position and
    force satisfy: eta_CO = H(r_bg)' * eta_CG and tau_CO = H(r_bg)' * tau_CG 
    """

    H = np.identity(6,float)
    H[0:3, 3:6] = Smtrx(r).T

    return H

L = 0.4571                # length (m)
#self.diam = 0.19            # cylinder diameter (m)
width = 0.436
height = 0.253

S = L * width
a = L/2                         # semi-axes
b = height/2                  
r_bg = np.array([0, 0, 0.02], float)    # CG w.r.t. to the CO
r_bb = np.array([0, 0, 0], float)       # CB w.r.t. to the CO

# Parasitic drag coefficient CD_0, i.e. zero lift and alpha = 0
# F_drag = 0.5 * rho * Cd * (pi * b^2)   
# F_drag = 0.5 * rho * CD_0 * S
Cd = 0.42                              # from Allen et al. (2000)
CD_0 = Cd * math.pi * b**2 / S

# Rigid-body mass matrix expressed in CO
m = 11.5     # mass of ROV
Ix = 0.16                      # moment of inertia
Iy = Ix
Iz = Ix
MRB_CG = np.diag([ m, m, m, Ix, Iy, Iz ])   # MRB expressed in the CG     
H_rg = Hmtrx(r_bg)
MRB = H_rg.T @ MRB_CG @ H_rg           # MRB expressed in the CO
g=9.8
# Weight and buoyancy
W = m * g
B = W

# Added moment of inertia in roll: A44 = r44 * Ix
r44 = 0.3           
MA_44 = r44 * Ix

# Lamb's k-factors (NOT SURE)
e = math.sqrt( 1-(b/a)**2 )
alpha_0 = ( 2 * (1-e**2)/pow(e,3) ) * ( 0.5 * math.log( (1+e)/(1-e) ) - e )  
beta_0  = 1/(e**2) - (1-e**2) / (2*pow(e,3)) * math.log( (1+e)/(1-e) )

k1 = alpha_0 / (2 - alpha_0)
k2 = beta_0  / (2 - beta_0)
k_prime = pow(e,4) * (beta_0-alpha_0) / ( 
    (2-e**2) * ( 2*e**2 - (2-e**2) * (beta_0-alpha_0) ) )   


# Added mass system matrix expressed in the CO
#MA = np.diag([ m*k1, m*k2, m*k2, MA_44, k_prime*Iy, k_prime*Iy ])
MA = np.diag([ m*k1, m*k2, m*k2, k_prime*Ix, k_prime*Iy, k_prime*Iy ])
print("MA matrix: \n", MA)