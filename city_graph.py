import matplotlib.pyplot as plt

def plot_city_graph(file_path):
    cities = []

    # Считываем данные из файла
    with open(file_path, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split())
            cities.append((x, y))

    if not cities:
        print("No cities to plot.")
        return

    # Разделяем координаты на X и Y
    x_coords = [city[0] for city in cities]
    y_coords = [city[1] for city in cities]

    # Создаём фигуру и оси
    plt.figure(figsize=(8, 6))

    # Рисуем города как точки
    plt.scatter(x_coords, y_coords, color="blue", label="Cities")

    # Рисуем путь между городами
    plt.plot(x_coords, y_coords, color="red", linestyle='-', marker='o', label="Tour Path")

    # Подписываем оси
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.title("City Graph - Tour Visualization")
    plt.legend()

    # Отображаем график
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    plot_city_graph("../.city_graph")
