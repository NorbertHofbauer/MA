import numpy as np
import matplotlib.pyplot as plt



# Read Data from file
temp_data1 = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/tf0_x0.res",
    skip_header=1
)
x = temp_data1[:, 0:1]  # 0
y = temp_data1[:, 1:2]  # 1
temp1 = temp_data1[:, 2:3]  # 2

# Read Data from file
temp_data2 = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/tf0_x1.res",
    skip_header=1
)
x = temp_data2[:, 0:1]  # 0
y = temp_data2[:, 1:2]  # 1
temp2 = temp_data2[:, 2:3]  # 2

error = np.zeros(len(y))
for i in range(len(y)):
 error[i]=(temp2[i] - temp1[i])
 
# Plot
fig, ax = plt.subplots(ncols=2,nrows=1,layout="constrained")
# Temperature field
ax[0].plot(y, temp1, '-')
ax[0].plot(y, temp2, '-.')
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
 temp_data[i,2]=temp1[i] #temp1
 temp_data[i,3]=temp2[i] #temp2
 temp_data[i,4]=error[i] #error
 temp_data[i,5]=error.mean() #error mean
 
np.savetxt("temp_advection.csv", temp_data, delimiter=",")

plt.show()
