import splinepy as spp
import numpy as np
import matplotlib.pyplot as plt

# Set Degree
degree = 4
#set refinement
ref = 4

# Recreate the geometry from file (as far as I understand)
geometry = spp.helpme.create.box(1, 1).bspline
geometry.elevate_degrees([0, 1] * degree)
geometry.insert_knots(0, np.linspace(0.125, 0.875, pow(2, ref)-1))
geometry.insert_knots(1, np.linspace(0.125, 0.875, pow(2, ref)-1))
# Note: this spline has the amazing property, that xi=x and y=eta, which means
# there is no need to do a "mapping" when performing the derivatives

# Vars
knot_vector = geometry.knot_vectors[0]
#knot_vector = [0, 0, 0, 0, 0.125, 0.125, 0.25, 0.25, 0.375, 0.375, 0.5, 0.5, 0.625, 0.625, 0.75, 0.75, 0.875, 0.875, 1, 1, 1, 1]

# Initiate Velocity Field
velocity_field = spp.BSpline(
    degrees=[degree],
    knot_vectors=[knot_vector],
    control_points=np.ones((len(knot_vector) - degree - 1, 2))
)

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

# Use splinepy to fit
# See documentation here:
# https://tataratat.github.io/splinepy/_generated/splinepy.helpme.fit.curve.html
_, residual = spp.helpme.fit.curve(
    fitting_points=uv,
    degree=degree,
    knot_vector=knot_vector,
    fitting_spline=velocity_field,
    associated_queries=y  # use y axis as parametric queries
)

print(f"Residual of the curve fit is {residual}")

# Evaluate first and second order derivatives with respect to parametric axis
u_v = velocity_field.evaluate(y)
du_dy = velocity_field.derivative(y, [1])
d2u_dy2 = velocity_field.derivative(y, [2])

# compute missing
K = 200000 * np.ones(len(kin_vis[:, 0]))
# compute d2u_dy2 with solver values
#print("d2u_dy2_solver")
#for i in range(len(kin_vis[:, 0])):
# print(-K[i],1/kin_vis[i,0].ravel())
#d2u_dy2_solver = np.multiply(-K[:],1/kin_vis[:,0].ravel())
d2u_dy2_solver = deriv[:,4]
# compute K_solver with solver values
print("K_solver")
for i in range(len(d2u_dy2_solver[:])):
 print(np.multiply(-d2u_dy2_solver[i].ravel(),kin_vis[i,0].ravel()))
K_solver = np.multiply(-d2u_dy2_solver[:].ravel(),kin_vis[:,0].ravel())


# compute dyn_vis_spp with splinepy values
dyn_vis_spp = np.zeros(len(y[:, 0])) #dummy

if visco_model[0] == 0:
 dyn_vis_spp = mp[0] * np.ones(len(y[:, 0]))
if visco_model[0] == 1:
 dyn_vis_spp = mp[0] / np.power(1 + mp[1] * abs(du_dy[:, 0]), mp[2])

print("dyn_vis_spp")
for i in range(len(dyn_vis_spp)):
 print(dyn_vis_spp[i])

# compute K with splinepy values
print("K_spp")
for i in range(len(d2u_dy2[:, 0])):
 print(np.multiply(-d2u_dy2[i].ravel(),dyn_vis_spp[i].ravel()))
K_spp = np.multiply(-d2u_dy2[:,0].ravel(),dyn_vis_spp[:].ravel())

#K_spp = np.gradient(dyn_vis_spp,y[:,0],edge_order=2)*du_dy[:,0] + dyn_vis_spp * d2u_dy2[:,0]
K_spp = np.gradient(dyn_vis_spp*du_dy[:,0],y[:,0],edge_order=2)
#K_spp = np.gradient(dyn_vis_spp,y[:,0],edge_order=2)*du_dy[:,0] + dyn_vis_spp * np.gradient(du_dy[:,0],y[:,0],edge_order=2)

# compute error 
error_u = np.multiply(uv[:, 0].ravel()-u_v[:, 0].ravel(),100/u_v[:, 0].ravel())
# check elementwise to catch zero values from u
print("u")
for i in range(len(u_v[:, 0])):
 print(u_v[i,0])
 if uv[i,0] == 0:
  error_u[i] = 0

error_du_dy = np.multiply(deriv[:, 1].ravel()-du_dy[:, 0].ravel(),100/du_dy[:, 0].ravel())
# check elementwise to catch zero values from du_dy
print("du_dy")
for i in range(len(du_dy[:, 0])):
 print(du_dy[i,0])
 if du_dy[i,0]*du_dy[i,0] < 1e-20:
  error_du_dy[i] = 0

error_d2u_dy2 = np.multiply(d2u_dy2_solver.ravel()-d2u_dy2[:, 0].ravel(),100/d2u_dy2[:, 0].ravel())
# check elementwise to catch zero values from du_dy
print("d2u_dy2")
for i in range(len(d2u_dy2[:, 0])):
 print(d2u_dy2[i,0])
 #if du_dy[i,0]*du_dy[i,0] < 1e-20:
 # error_du_dy[i] = 0

error_eta = np.multiply(kin_vis[:,0].ravel()-dyn_vis_spp[:].ravel(),100/dyn_vis_spp[:].ravel())
print("eta")
for i in range(len(dyn_vis_spp)):
 print(dyn_vis_spp[i])

error_K = np.multiply(K_solver-K_spp.ravel(),100/K_spp.ravel())
# check elementwise to catch zero values from du_dy
print("-K")
for i in range(len(K_spp)):
 print(-K_spp[i])
 #if du_dy[i,0]*du_dy[i,0] < 1e-20:
 # error_du_dy[i] = 0


# Plot
fig, ax = plt.subplots(ncols=3,nrows=5,layout="constrained")
# Velocity field
ax[0,0].plot(y, u_v[:, 0].ravel(), '-x')
ax[0,0].set_title("Velocity field properties")
ax[0,0].set_ylabel("Velocity")
ax[0,0].set_xlabel("y-axis")
# First derivative field
ax[1,0].plot(y, du_dy[:, 0].ravel(), '-x')
ax[1,0].set_ylabel("Shear Rate\n(du/dy)")
ax[1,0].set_xlabel("y-axis")
# Second derivative field
ax[2,0].plot(y, d2u_dy2[:, 0].ravel(), '-x')
ax[2,0].set_ylabel("Stress-factor\n(d2u/dy2)")
ax[2,0].set_xlabel("y-axis")
# compute eta
ax[3,0].plot(y, dyn_vis_spp[:].ravel(), '-x')
ax[3,0].set_ylabel("eta")
ax[3,0].set_xlabel("y-axis")
# compute eta*d2u_dy2
#plt.plot(y, np.multiply(d2u_dy2[:, 0].ravel(),kin_vis), '-x')
ax[4,0].plot(y, -K_spp[:].ravel(), '-x')
ax[4,0].set_ylabel("eta * (d2u/dy2) = -K")
ax[4,0].set_xlabel("y-axis")

# Velocity field
ax[0,1].plot(y, uv[:, 0].ravel(), '-x')
ax[0,1].set_ylabel("Velocity")
ax[0,1].set_xlabel("y-axis")
# First derivative field
ax[0,1].set_title("solver")
ax[1,1].plot(y, deriv[:, 1].ravel(), '-x')
ax[1,1].set_ylabel("Shear Rate\n(du/dy)")
ax[1,1].set_xlabel("y-axis")
# Second derivative field
ax[2,1].plot(y, d2u_dy2_solver[:].ravel(), '-x')
ax[2,1].set_ylabel("Stress-factor \n(d2u/dy2)")
ax[2,1].set_xlabel("y-axis")
# compute eta
ax[3,1].plot(y, kin_vis[:,0], '-x')
ax[3,1].set_ylabel("eta")
ax[3,1].set_xlabel("y-axis")
# compute eta*d2u_dy2
ax[4,1].plot(y, -K_solver, '-x')
ax[4,1].set_ylabel("eta * (d2u/dy2) = -K")
ax[4,1].set_xlabel("y-axis")

# Velocity field
ax[0,2].set_title("error")
ax[0,2].plot(y, error_u, '-x')
ax[0,2].set_ylabel("Velocity u\nrel error %")
ax[0,2].set_xlabel("y-axis")
# First derivative field
ax[1,2].plot(y, error_du_dy, '-x')
ax[1,2].set_ylabel("Shear Rate\nrel error %")
ax[1,2].set_xlabel("y-axis")
# Second derivative field
ax[2,2].plot(y, error_d2u_dy2, '-x')
ax[2,2].set_ylabel("d2u_dy2\nrel error %")
ax[2,2].set_xlabel("y-axis")
# compute eta
ax[3,2].plot(y, error_eta[:].ravel(), '-x')
ax[3,2].set_ylabel("eta\nrel error %")
ax[3,2].set_xlabel("y-axis")
# compute eta*d2u_dy2
ax[4,2].plot(y, error_K[:].ravel(), '-x')
ax[4,2].set_ylabel("eta * (d2u/dy2)\nrel error %")
ax[4,2].set_xlabel("y-axis")

plt.show()

#print("0 0 0 0 0.125 0.125 0.25 0.25 0.375 0.375 0.5 0.5 0.625 0.625 0.75 0.75 0.875 0.875 1 1 1 1")
#print(knot_vector)
