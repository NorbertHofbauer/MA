import sys, getopt
import numpy
import matplotlib.pyplot as plt

def read_scalar(filepath):
    with open(filepath) as fp:
        line = fp.readline()
        name = line
        line = fp.readline()
        cnt = 1
        data_array = []
        while line:
            print("Line {}: {}".format(cnt, line.strip()))
            data_array.append([float(x) for x in line.split()])
            line = fp.readline()
            cnt += 1
        #print(name, " scalar ",data_array)
        plot_scalar(data_array,name)

def read_vector(filepath):
    with open(filepath) as fp:
        line = fp.readline()
        name = line
        line = fp.readline()
        cnt = 1
        data_array = []
        while line:
            print("Line {}: {}".format(cnt, line.strip()))
            data_array.append([float(x) for x in line.split()])
            line = fp.readline()
            cnt += 1
        #print(name, "vector ",data_array)
        plot_vector(data_array,name)

def read_fsi(filepath):
    with open(filepath) as fp:
        line = fp.readline()
        name = line
        line = fp.readline()
        cnt = 1
        data_array = []
        while line:
            print("Line {}: {}".format(cnt, line.strip()))
            data_array.append([float(x) for x in line.split()])
            line = fp.readline()
            cnt += 1
        #print(name, " fsi scalar ",data_array)
        plot_scalar_fsi(data_array,name)

def plot_scalar(data, name_x_axis):
    x_data = []
    y_data = []
    for item in data:
        x_data.append(item[2])
        y_data.append(item[1])
    x = numpy.array(x_data)
    y = numpy.array(y_data)
    plt.plot(x,y, marker="o")
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.ylabel('y')
    plt.xlabel(name_x_axis)
    plt.show()

def plot_vector(data, name_x_axis):
    x_data = []
    y_data = []
    for item in data:
        x_data.append(item[2])
        y_data.append(item[1])
    x = numpy.array(x_data)
    y = numpy.array(y_data)
    plt.plot(x,y, marker="o")
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.ylabel('y')
    plt.xlabel("x " + name_x_axis)
    plt.show()

    x_data = []
    y_data = []
    for item in data:
        x_data.append(item[3])
        y_data.append(item[1])
    x = numpy.array(x_data)
    y = numpy.array(y_data)
    plt.plot(x,y, marker="o")
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.ylabel('y')
    plt.xlabel("y " + name_x_axis)
    plt.show()

def plot_scalar_fsi(data, name_y_axis):
    x_data = []
    y_data = []
    for item in data:
        x_data.append(item[0])
        y_data.append(item[2])
    x = numpy.array(x_data)
    y = numpy.array(y_data)
    plt.plot(x,y, marker="o")
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.xlabel('x')
    plt.ylabel(name_y_axis)
    plt.show()

def main(argv):
    input_folder = ''
    opts, args = getopt.getopt(argv,"hi:",["inputfolder="])
    for opt, arg in opts:
        if opt == '-h':
            print ('postprocess.py -i <inputfolder>')
            sys.exit()
        elif opt in ("-i", "--inputfolder"):
            input_folder = arg
    print ('Input tf file is ', input_folder + "tf.res")
    print ('Input v file is ', input_folder + "v.res")
    print ('Input fsi file is ', input_folder + "fsi.res")

    if not input_folder == '':
        read_scalar(input_folder + "tf.res")
        read_vector(input_folder + "v.res")
        read_fsi(input_folder + "fsi.res")

if __name__ == "__main__":
    main(sys.argv[1:])