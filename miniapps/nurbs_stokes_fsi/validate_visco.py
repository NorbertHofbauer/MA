import numpy as np
import matplotlib.pyplot as plt


# Read Data from file
visco_model_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/visco.res",
    skip_header=1,max_rows=1
)
visco_model = visco_model_data[0:1]
mp = visco_model_data[1:6]


# Read Data from file
file_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/visco.res",
    skip_header=3  # First three lines are supplementary info
)
x = file_data[:, 0:1]  # 0
y = file_data[:, 1:2]  # 1
uv = file_data[:, 2:4]  # 2,3
kin_vis = file_data[:, 4:5]  # 4
gen_shear = file_data[:, 5:6]  # 5
deriv = file_data[:, 6:11]  # 6,7,8,9,10 all derivates 

# Read Data from file
temp_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/tf0.res",
    skip_header=1
)
temp_x = temp_data[:, 0:1]  # 0
temp_y = temp_data[:, 1:2]  # 1
temp = temp_data[:, 2:3]  # 2

#d2u_dy2_solver = np.multiply(-K[:],1/kin_vis[:,0].ravel())
d2u_dy2_solver = deriv[:,4]
# compute K_solver with solver values
K_solver = np.gradient(kin_vis[:,0] * deriv[:,1],y[:,0],edge_order=2)
print("K_solver")
for i in range(len(K_solver[:])):
 print(K_solver[i])


#error_K = np.multiply(K_solver-K_spp[:,0].ravel(),100/K_spp[:,0].ravel())
# check elementwise to catch zero values from du_dy
#print("-K")
#for i in range(len(K_spp)):
# print(-K_spp[i])
 #if du_dy[i,0]*du_dy[i,0] < 1e-20:
 # error_du_dy[i] = 0


# Plot
#fig, ax = plt.subplots(ncols=3,nrows=5,layout="constrained")
fig, ax = plt.subplots(ncols=2,nrows=3,layout="constrained")
# Velocity field
ax[0,0].plot(y, uv[:, 0].ravel(), '-')
ax[0,0].set_title("Velocity field properties")
ax[0,0].set_ylabel("Velocity")
ax[0,0].set_xlabel("y-axis")
# First derivative field
ax[1,0].plot(y, deriv[:, 1].ravel(), '-')
ax[1,0].set_ylabel("Shear Rate\n(du/dy)")
ax[1,1].set_xlabel("y-axis")
# Second derivative field
ax[2,0].plot(y, d2u_dy2_solver[:].ravel(), '-')
ax[2,0].set_ylabel("Stress-factor\n(d2u/dy2)")
ax[2,0].set_xlabel("y-axis")
# temperature field
ax[0,1].plot(y, temp[:,0].ravel(), '-')
ax[0,1].set_ylabel("Temperature")
ax[0,1].set_xlabel("y-axis")
# eta
ax[1,1].plot(y, kin_vis[:,0].ravel(), '-')
ax[1,1].set_ylabel("eta")
ax[1,1].set_xlabel("y-axis")
# d(eta*d2u_dy2)/dy
ax[2,1].plot(y, K_solver[:].ravel(), '-x')
ax[2,1].plot(y, np.ones(len(K_solver))*K_solver[:].mean(), '-')
print('K mean value')
print(K_solver[:].mean())
print('K norm')
print(np.linalg.norm(K_solver[:],np.inf))
ax[2,1].set_ylabel("-K")
ax[2,1].set_title("mean: " + str(K_solver[:].mean()))
ax[2,1].set_xlabel("y-axis")


plt.show()

