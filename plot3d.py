import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re

def read_data_from_plot3d_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    x_list = []
    y_data = []
    _x_list = []
    _y_data = []

    # Парсим данные
    for line in lines:
        if line.startswith("x"):
            index = int(re.search(r"x(\d+):", line).group(1))
            if len(x_list) <= index:
                x_list.append([])
            x_values = re.findall(r'-?\d+\.\d+|-?\d+', line.split(":")[1])
            x_list[index] = [float(value) for value in x_values]
        elif line.startswith("y:"):
            y_values = re.findall(r'-?\d+\.\d+|-?\d+', line[2:])
            y_data = [float(value) for value in y_values]
        elif line.startswith("_x"):
            index = int(re.search(r"_x(\d+):", line).group(1))
            if len(_x_list) <= index:
                _x_list.append([])
            x_values = re.findall(r'-?\d+\.\d+|-?\d+', line.split(":")[1])
            _x_list[index] = [float(value) for value in x_values]
        elif line.startswith("_y:"):
            y_values = re.findall(r'-?\d+\.\d+|-?\d+', line[3:])
            _y_data = [float(value) for value in y_values]

    return x_list, y_data, _x_list, _y_data

def plot_3d_and_2d_projections(x_list, y_data, _x_list, _y_data):
    fig = plt.figure(figsize=(16, 10))

    if len(x_list) == 2:
        X, Y = np.meshgrid(np.unique(x_list[0]), np.unique(x_list[1]))
        Z = np.array(y_data).reshape(X.shape)

        # 3D график
        ax_3d = fig.add_subplot(221, projection='3d')
        ax_3d.plot_surface(X, Y, Z, cmap='viridis', edgecolor='k', alpha=0.8)
        if len(_x_list) == 2 and len(_y_data) > 0:
            ax_3d.scatter(_x_list[0], _x_list[1], _y_data, label="Population", color='red', zorder=5)

        ax_3d.set_xlabel("X1")
        ax_3d.set_ylabel("X2")
        ax_3d.set_zlabel("Y")
        ax_3d.set_title("3D Plot")
        ax_3d.legend()

        # 2D проекции
        # Проекция на X1-Y
        ax_xy = fig.add_subplot(222)
        ax_xy.plot(x_list[0], y_data, label="Projection X1-Y", color='green')
        if len(_x_list) == 2 and len(_y_data) > 0:
            ax_xy.scatter(_x_list[0], _y_data, label="Population X1-Y", color='red')
        ax_xy.set_xlabel("X1")
        ax_xy.set_ylabel("Y")
        ax_xy.set_title("Projection X1-Y")
        ax_xy.legend()

        # Проекция на X2-Y
        ax_xz = fig.add_subplot(223)
        ax_xz.plot(x_list[1], y_data, label="Projection X2-Y", color='purple')
        if len(_x_list) == 2 and len(_y_data) > 0:
            ax_xz.scatter(_x_list[1], _y_data, label="Population X2-Y", color='red')
        ax_xz.set_xlabel("X2")
        ax_xz.set_ylabel("Y")
        ax_xz.set_title("Projection X2-Y")
        ax_xz.legend()

        # Проекция на X1-X2
        ax_xx = fig.add_subplot(224)
        ax_xx.plot(x_list[0], x_list[1], label="Projection X1-X2", color='orange')
        if len(_x_list) == 2 and len(_y_data) > 0:
            ax_xx.scatter(_x_list[0], _x_list[1], label="Population X1-X2", color='red')
        ax_xx.set_xlabel("X1")
        ax_xx.set_ylabel("X2")
        ax_xx.set_title("Projection X1-X2")
        ax_xx.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    file_path = "../.plot3d"
    x_list, y_data, _x_list, _y_data = read_data_from_plot3d_file(file_path)
    plot_3d_and_2d_projections(x_list, y_data, _x_list, _y_data)
