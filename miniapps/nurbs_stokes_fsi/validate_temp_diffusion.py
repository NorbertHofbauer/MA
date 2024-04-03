import numpy as np
import matplotlib.pyplot as plt



# Read Data from file
file_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/tf0.res",
    skip_header=1  # First line 
)
x = file_data[:, 0:1]  # 0
y = file_data[:, 1:2]  # 1
temp = file_data[:, 2:3]  # 2

# Read Data from file
temp_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/tf0.res",
    skip_header=1
)
x = temp_data[:, 0:1]  # 0
y = temp_data[:, 1:2]  # 1
temp = temp_data[:, 2:3]  # 2

Tmin=100
Tmax=200

analytic = np.linspace(Tmin, Tmax, len(y), endpoint=True)

error = np.zeros(len(y))
for i in range(len(y)):
 error[i]=(temp[i] - analytic[i])
 
# Plot
fig, ax = plt.subplots(ncols=2,nrows=1,layout="constrained")
# Temperature field
ax[0].plot(y, temp, '-')
ax[0].plot(y, analytic, '-.')
ax[0].set_title("Temperature field")
ax[0].set_ylabel("Temperature")
ax[0].set_xlabel("y-axis")
# Error
ax[1].plot(y, error, '-')
ax[1].set_title("Absolute Error")
#ax[1].set_ylabel("absolute Error")
ax[1].set_xlabel("y-axis")

temp_data = np.zeros((len(y),6))
for i in range(len(y)):
 temp_data[i,0]=x[i] #x coord
 temp_data[i,1]=y[i] #y coord
 temp_data[i,2]=temp[i] #temp
 temp_data[i,3]=analytic[i] #analytic
 temp_data[i,4]=error[i] #error
 temp_data[i,5]=error.mean() #error mean
 
np.savetxt("temp_diffusion.csv", temp_data, delimiter=",")

plt.show()
