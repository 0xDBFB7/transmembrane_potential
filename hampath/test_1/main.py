# Main file for the simple shooting case
#
# Problem definition
# min int_0^1 u^2
# dot(x) = v
# dox(v) = -lambda v^2 + u
# x(0) = x_0, x(1) = x_f, v(0) = v_0, v(1) = v_f

import numpy as np
import matplotlib.pyplot as plt

# It is assumed that HamPath lib files are placed in libhampath
from libhampath.hampathcode import *   # Change directory name if needed
from libhampath.control_p import *     # Change directory name if needed

#-------------------------------------------------------------------------------------------------------------%
print('\nStep 1: parameters initialization\n')

t0      = 0.0                              # Initial time
tf      = 1.0                              # Final time
q0      = np.array([-1.0,0.0])             # Initial state
n       = len(q0)                          # State dimension
yGuess  = np.array([0.1,0.1])              # Initial guess for the shooting metdhod: yGuess = [px(0) pv(0)]
par     = np.array([t0, tf, q0[0], q0[1], 0.0, 0.0, 0.0])  # t0, tf, x_0, v_0, x_f, v_f, lambda_0
npar    = len(par)
par0    = par                              # par_0
parf    = np.array(par); parf[-1] = 1.0    # par_f with lambda_f = 1.0: for homotopy on lambda
options = HampathOptions()                 # Hampath options
print(options)

#-------------------------------------------------------------------------------------------------------------%
print('\nStep 2: first shooting\n')

[y0,ssol,nfev,njev,flag] = ssolve(yGuess,options,par0)

#-------------------------------------------------------------------------------------------------------------%
print('\nStep 3: continuation on lambda from ', par0[-1], ' to ', parf[-1], '\n')

parspan = np.zeros((npar,2))
parspan[:,0] = par0
parspan[:,1] = parf

[ parout, yout, sout, viout, dets, normS, ps, flag ] = hampath(parspan,y0,options)

lout = parout[-1,:]
p0f  = yout[:,-1]
parf = parout[:,-1]

# Figures
fig, axarr = plt.subplots(nrows=2, ncols=4)
#fig.tight_layout() # Or equivalently,  "plt.tight_layout()"

left    = 0.125 # the left side of the subplots of the figure
right   = 0.9   # the right side of the subplots of the figure
bottom  = 0.1   # the bottom of the subplots of the figure
top     = 0.9   # the top of the subplots of the figure
wspace  = 0.3   # the amount of width reserved for blank space between subplots
hspace  = 0.3   # the amount of height reserved for white space between subplots

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

# Paths
lig,col = 0,0
axarr[lig,col].plot(yout[0,:],lout,'k')
axarr[lig,col].set_xlabel('$p_x(0)$')
axarr[lig,col].set_ylabel('$\lambda$')
axarr[lig,col].set_title('Paths')

lig,col = 1,0
axarr[lig,col].plot(yout[1,:],lout,'k')
axarr[lig,col].set_xlabel('$p_v(0)$')
axarr[lig,col].set_ylabel('$\lambda$')

# Solution
z0          = np.zeros((2*n,))
z0[0:n]     = q0
z0[n:2*n]   = p0f
[ tout, z, flag ] = exphvfun(np.array([t0, tf]), z0, options, parf)

lig,col = 0,1
axarr[lig,col].plot(tout,z[0,:],'b')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_ylabel('$x(t)$')
axarr[lig,col].set_title('State solution')

lig,col = 0,2
axarr[lig,col].plot(tout,z[2,:],'b')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_ylabel('$p_x(t)$')
axarr[lig,col].set_title('Co-state solution')

lig,col = 1,1
axarr[lig,col].plot(tout,z[1,:],'b')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_ylabel('$v(t)$')

lig,col = 1,2
axarr[lig,col].plot(tout,z[3,:],'b')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_ylabel('$p_v(t)$')

# Control
u   = control(tout,z,parf)
lig,col = 0,3
axarr[lig,col].plot(tout,u[0,:],'r')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_ylabel('$u(t)$')
axarr[lig,col].set_title('Control')

#-------------------------------------------------------------------------------------------------------------%
print('\nStep 4: optimality checked on final solution\n')

k = n # regular case, k = n : fixed tf
dz0 = np.zeros((2*n,n))
dz0[n:2*n,:] = np.eye(n)

[ tout, z, dz, flag ] = expdhvfun(np.array([t0, tf]), z0, dz0, options, parf)

nt  = len(tout)
sv  = np.zeros(nt,)
de  = np.zeros(nt,)

for j in range(0,nt):
    dq      = dz[0:n,j*k:(j+1)*k]
    sv[j]   = np.min(np.linalg.svd(dq,compute_uv=0))    # Get smallest singular value
    de[j]   = np.linalg.det(dq)                         # Get determinant

maxsv       = np.max(sv)
lig,col     = 1,3
plt1,       = axarr[lig,col].plot(tout,sv,'m',label='$\sigma_{min}$')
plt2,       = axarr[lig,col].plot(tout,maxsv*np.sign(de),'k--',label='$sign(\det(\delta q))$')
axarr[lig,col].set_xlabel('t')
axarr[lig,col].set_title('Conjugate points')
plt.legend([plt1, plt2])

#-------------------------------------------------------------------------------------------------------------%
plt.show()
