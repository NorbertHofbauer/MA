import splinepy as spp
import numpy as np
import matplotlib.pyplot as plt

# Set Degree
degree = 3

# Recreate the geometry from file (as far as I understand)
geometry = spp.helpme.create.box(1, 1).bspline
geometry.elevate_degrees([0, 1] * degree)
geometry.insert_knots(0, np.linspace(0.125, 0.875, 7))
geometry.insert_knots(1, np.linspace(0.125, 0.875, 7))
# Note: this spline has the amazing property, that xi=x and y=eta, which means
# there is no need to do a "mapping" when performing the derivatives

# Vars
knot_vector = geometry.knot_vectors[0]

# Initiate Velocity Field
velocity_field = spp.BSpline(
    degrees=[degree],
    knot_vectors=[knot_vector],
    control_points=np.ones((len(knot_vector) - degree - 1, 2))
)


# Read Data from file
file_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/visco.res",
    skip_header=2  # First two lines are supplementary info
)
x = file_data[:, 0:1]  # 0
y = file_data[:, 1:2]  # 1
uv = file_data[:, 2:4]  # 2,3
kin_vis = file_data[:, 4:5]  # 4
gen_shear = file_data[:, 5:6]  # 5
deriv = file_data[:, 6:10]  # 6,7,8,9 all derivates 

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
du_dy = velocity_field.derivative(y, [1])
d2u_dy2 = velocity_field.derivative(y, [2])

# compute error 
error = np.multiply(deriv[:, 1].ravel()-du_dy[:, 0].ravel(),100/du_dy[:, 0].ravel())
# check elementwise to catch zero values from du_dy
for i in range(len(du_dy[:, 0])):
 print(du_dy[i,0])
 if du_dy[i,0]*du_dy[i,0] < 1e-20:
  error[i] = 0

for i in range(len(d2u_dy2[:, 0])):
 print(np.multiply(d2u_dy2[i, 0].ravel(),kin_vis[i,0].ravel()))


# Plot
fig, ax = plt.subplots(ncols=3,nrows=4,layout="constrained")
# Velocity field
ax[0,0].plot(y, uv[:, 0].ravel(), '-x')
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
# compute eta*d2u_dy2
#plt.plot(y, np.multiply(d2u_dy2[:, 0].ravel(),kin_vis), '-x')
ax[3,0].plot(y, np.multiply(d2u_dy2[:, 0].ravel(),kin_vis[:,0].ravel()), '-x')
ax[3,0].set_ylabel("eta * (d2u/dy2)")
ax[3,0].set_xlabel("y-axis")

# First derivative field
ax[0,1].set_title("solver")
ax[1,1].plot(y, deriv[:, 1].ravel(), '-x')
ax[1,1].set_ylabel("Shear Rate\n(du/dy)")
ax[1,1].set_xlabel("y-axis")

# First derivative field
ax[0,2].set_title("error")
ax[1,2].plot(y, error, '-x')
ax[1,2].set_ylabel("Shear Rate\nrel error %")
ax[1,2].set_xlabel("y-axis")

plt.show()

