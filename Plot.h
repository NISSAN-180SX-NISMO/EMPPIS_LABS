#ifndef PLOT_H
#define PLOT_H

#include "vector"
#include "fstream"

inline std::ofstream out;

namespace plot
{

template<typename T1, typename T2, typename T3, typename T4>
inline void print(std::vector<T1> x, std::vector<T2> y, std::vector<T3> _x = {}, std::vector<T4> _y = {})
{
    // x, y  оригинальная функция, где х это обычно пчисла равномерно растянутые на промежутке какомто
    // а у это значения некой функции в этих точках

    // _x и _y это моя текйщая популяция, просто какието значения (особи), и значение функции в этих точках
    out.open("../.plot");
    std::string _out;
    _out = "x:";
    for (const auto& _ : x)
        _out += std::to_string(float(_)) + ",";

    _out += "\ny:";
    for (const auto& _ : y)
        _out += std::to_string(float(_)) + ",";

    _out += "\n_x:";
    for (const auto& _ : _x)
        _out += std::to_string(float(_)) + ",";

    _out += "\n_y:";
    for (const auto& _ : _y)
        _out += std::to_string(float(_)) + ",";

    out << _out;
    out.close();
    system("C:/Users/User/Desktop/python/SGA/.venv/Scripts/python.exe ../plot.py");
}


template<typename T1, typename T2, typename T3, typename T4>
inline void print3d(std::vector<std::vector<T1>> x_list, std::vector<T2> y, std::vector<std::vector<T3>> _x_list = {}, std::vector<T4> _y = {})
{
    // Открываем файл для записи
    out.open("../.plot3d");
    std::string _out;

    // Записываем x_list (каждый x как отдельный набор значений)
    for (size_t i = 0; i < x_list.size(); ++i)
    {
        _out += "x" + std::to_string(i) + ":";
        for (const auto& val : x_list[i])
        {
            _out += std::to_string(float(val)) + ",";
        }
        _out += "\n";
    }

    // Записываем y
    _out += "y:";
    for (const auto& val : y)
    {
        _out += std::to_string(float(val)) + ",";
    }
    _out += "\n";

    // Записываем _x_list
    std::vector<std::string> _out_list;
    if (not _x_list.empty())
        _out_list.resize(_x_list.front().size());

    for (const auto& pair: _x_list)
    {
        for (size_t i = 0; i < pair.size(); ++i)
            _out_list[i] += std::to_string(float(pair[i])) + ",";
    }

    for (size_t i = 0; i < _out_list.size(); ++i)
    {
        _out += "_x" + std::to_string(i) + ":";
        for (const auto& val : _out_list[i])
        {
            _out += val;
        }
        _out += "\n";
    }

    // Записываем _y
    _out += "_y:";
    for (const auto& val : _y)
    {
        _out += std::to_string(float(val)) + ",";
    }
//    _out.erase(_out.size() - 1);
    _out += "\n";

    out << _out;
    out.close();

    // Вызываем Python-скрипт для построения графика
    system("C:/Users/User/Desktop/python/SGA/.venv/Scripts/python.exe ../plot3d.py");
}

inline void city_graph(const std::vector<std::pair<float, float>> &tour)
{
    std::ofstream out("../.city_graph");

    if (!out.is_open())
    {
        throw std::runtime_error("Failed to open .city_graph for writing.");
    }

    // Записываем города в файл. Каждая пара координат с новой строки
    for (const auto &city : tour)
    {
        out << city.first << " " << city.second << "\n";
    }

    out.close();

    // Вызываем Python-скрипт для построения графика
    system("C:/Users/User/Desktop/python/SGA/.venv/Scripts/python.exe ../city_graph.py");
}

}

#endif //PLOT_H
