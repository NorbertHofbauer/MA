import numpy as np
import matplotlib.pyplot as plt



# Read Data from file
file_data = np.genfromtxt(
    "/home/user/MA/build/miniapps/nurbs_stokes_fsi/fsi.res",
    skip_header=1
)
x = file_data[:, 0:1]  # 0
y = file_data[:, 1:2]  # 1
temp = file_data[:, 2:3]  # 2

# split into parts
nov = int(len(temp)/3) # number of values

x1 = x[0:nov]
y1 = y[0:nov]
temp1 = temp[0:nov]
x2 = x[nov:2*nov]
y2 = y[nov:2*nov]
temp2 = temp[nov:2*nov]
x3 = x[2*nov:3*nov]
y3 = y[2*nov:3*nov]
temp3 = temp[2*nov:3*nov]


#error = np.zeros(len(y))
#for i in range(len(y)):
# error[i]=(temp2[i] - temp1[i])
 
# Plot
fig, ax = plt.subplots(ncols=2,nrows=1,layout="constrained")
# Temperature field
ax[0].plot(x1, temp1, '-')
ax[0].plot(x2, temp2, '-')
ax[0].plot(x3, temp3, '-')
ax[0].set_title("Temperature field")
ax[0].set_ylabel("Temperature")
ax[0].set_xlabel("x-axis")
# Error
#ax[1].plot(y, error, '-')
#ax[1].set_title("Absolute Error")
#ax[1].set_ylabel("absolute Error")
#ax[1].set_xlabel("y-axis")

temp_data = np.zeros((nov,6))
for i in range(nov):
 temp_data[i,0]=x1[i] #x1 coord
 temp_data[i,1]=temp1[i] #temp1
 temp_data[i,2]=x2[i] #x2 coord
 temp_data[i,3]=temp2[i] #temp2
 temp_data[i,4]=x3[i] #x3 coord
 temp_data[i,5]=temp3[i] #temp3
 
np.savetxt("temp_fsi.csv", temp_data, delimiter=",")

plt.show()
