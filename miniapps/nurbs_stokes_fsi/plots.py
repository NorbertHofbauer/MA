#!python

import matplotlib.pyplot as plt
import numpy as np

y = np.array([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])

x1 = np.array([0,0.975,1.9,2.775,3.6,4.375,5.1,5.775,6.4,6.975,7.5,7.975,8.4,8.775,9.1,9.375,9.6,9.775,9.9,9.975,10])
x2 = np.array([0,0.7375,1.45,2.1375,2.8,3.4375,4.05,4.6375,5.2,5.7375,6.25,6.7375,7.2,7.6375,8.05,8.4375,8.8,9.1375,9.45,9.7375,10])
x3 = np.array([0,0.025,0.1,0.225,0.4,0.625,0.9,1.225,1.6,2.025,2.5,3.025,3.6,4.225,4.9,5.625,6.4,7.225,8.1,9.025,10])
x4 = np.array([0,0.2625,0.55,0.8625,1.2,1.5625,1.95,2.3625,2.8,3.2625,3.75,4.2625,4.8,5.3625,5.95,6.5625,7.2,7.8625,8.55,9.2625,10])
x5 = np.array([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10])
x6 = np.array([0,0.475,0.9,1.275,1.6,1.875,2.1,2.275,2.4,2.475,2.5,2.475,2.4,2.275,2.1,1.875,1.6,1.275,0.9,0.475,0])


fig, ax = plt.subplots(layout="constrained")

ax.plot(x1, y, 'x--',color='b',markersize=8, label=r'K=20, U=10, $\eta$=1')
ax.plot(x1, y, 'o-',fillstyle='none',color='b',alpha=0.5, label=r'K=20, U=10, $\eta$=1, Theorie')
ax.plot(x2, y, 'x--',color='g',markersize=8, label=r'K=20, U=10, $\eta$=2')
ax.plot(x2, y, 'o-',fillstyle='none',color='g',alpha=0.5, label=r'K=20, U=10, $\eta$=2, Theorie')
ax.plot(x3, y, 'x--',color='c',markersize=8, label=r'K=-20, U=10, $\eta$=1')
ax.plot(x3, y, 'o-',fillstyle='none',color='c',alpha=0.5, label=r'K=-20, U=10, $\eta$=1, Theorie')
ax.plot(x4, y, 'x--',color='m',markersize=8, label=r'K=-20, U=10, $\eta$=2')
ax.plot(x4, y, 'o-',fillstyle='none',color='m',alpha=0.5, label=r'K=-20, U=10, $\eta$=2, Theorie')
ax.plot(x5, y, 'x--',color='gold',markersize=8, label=r'K=0, U=10, $\eta$=1')
ax.plot(x5, y, 'o-',fillstyle='none',color='gold',alpha=0.5, label=r'K=0, U=10, $\eta$=1, Theorie')
ax.plot(x6, y, 'x--',color='olive',markersize=8, label=r'K=20, U=0, $\eta$=1')
ax.plot(x6, y, 'o-',fillstyle='none',color='olive',alpha=0.5, label=r'K=20, U=0, $\eta$=1, Theorie')
ax.legend()
ax.legend(bbox_to_anchor=(0., -0.02, 1., -0.102), loc='upper left', ncols=2, mode="expand", borderaxespad=0.)
ax.set_title('Couette-Poiseuille')
ax.set_xlabel('u(y)')
ax.set_ylabel('y')
plt.grid()
plt.show()
