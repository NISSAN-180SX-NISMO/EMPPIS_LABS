import matplotlib.pyplot as plt
import re

def read_data_from_plot_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    x_data = []
    y_data = []
    _x_data = []
    _y_data = []


    # Parsing x and y data
    for line in lines:
        if line.startswith("x:"):
            x_values = re.findall(r'-?\d+\.\d+|-?\d+', line[2:])
            x_data = [float(value) for value in x_values]
        elif line.startswith("y:"):
            y_values = re.findall(r'-?\d+\.\d+|-?\d+', line[2:])
            y_data = [float(value) for value in y_values]
        elif line.startswith("_x:"):
            x_values = re.findall(r'-?\d+\.\d+|-?\d+', line[2:])
            _x_data = [float(value) for value in x_values]
        elif line.startswith("_y:"):
            y_values = re.findall(r'-?\d+\.\d+|-?\d+', line[2:])
            _y_data = [float(value) for value in y_values]

    return x_data, y_data, _x_data, _y_data

def plot_data(x_data, y_data, _x_data, _y_data):
    plt.plot(x_data, y_data, color='blue')
    plt.scatter(_x_data, _y_data, color='red', zorder=5)
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.title('Plot of X vs Y')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    file_path = "C:/Users/User/Desktop/cpp/gen_alg_lab1/.plot"
    x_data, y_data, _x_data, _y_data = read_data_from_plot_file(file_path)
    plot_data(x_data, y_data, _x_data, _y_data)