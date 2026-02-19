import numpy as np                                        
data = np.loadtxt('output/default00_cl.dat')
ell = data[:, 0]
Cl_TT = data[:, 1]
Dl_TT = Cl_TT * (2.7255e6)**2

for i in range(9):
    print(f"l={int(ell[i]):3d}  D_l={Dl_TT[i]:.2f} muK^2")
