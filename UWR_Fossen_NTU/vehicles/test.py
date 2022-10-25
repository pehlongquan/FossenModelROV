import numpy as np
import math
import sys





L = 0.4571                # length (m)
#self.diam = 0.19            # cylinder diameter (m)
width = 0.436
height = 0.253

a = L/2                         # semi-axes
b = height/2    

e = math.sqrt( 1-(b/a)**2 )
alpha_0 = ( 2 * (1-e**2)/pow(e,3) ) * ( 0.5 * math.log( (1+e)/(1-e) ) - e )  
beta_0  = 1/(e**2) - (1-e**2) / (2*pow(e,3)) * math.log( (1+e)/(1-e) )

k1 = alpha_0 / (2 - alpha_0)
k2 = beta_0  / (2 - beta_0)
k_prime = pow(e,4) * (beta_0-alpha_0) / ( 
            (2-e**2) * ( 2*e**2 - (2-e**2) * (beta_0-alpha_0) ) ) 


print(k1)
print(k2)
print(k_prime)
        
        