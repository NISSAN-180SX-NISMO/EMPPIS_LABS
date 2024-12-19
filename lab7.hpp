#ifndef LAB7_HPP
#define LAB7_HPP


#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include "Bocchiemark.hpp"

namespace lab7
{
class Swarm;
class Bug;
// Пример функции для минимизации
namespace k
{
uint16_t N = 100;
const float MIN = -5.12, MAX = 5.12;
const float VEL_MIN = -1.0; // Минимальная скорость
const float VEL_MAX = 1.0;  // Максимальная скорость
const float C1 = 2.0;       // Коэффициент индивидуального ускорения
const float C2 = 2.0;       // Коэффициент социального ускорения
const float INERTIA = 0.5;  // Инерционная компонента
const size_t MAX_ITERATIONS = 100; // Максимальное количество итераций

inline float f(const std::vector<float> &x);
inline std::vector<float> f(const Swarm &x);
}  // namespace k

// Класс частицы
class Bug {
public:
    std::vector<float> position;     // Текущая позиция
    std::vector<float> velocity;     // Текущая скорость
    std::vector<float> best_position; // Лучшая найденная позиция
    float best_value;                // Значение функции в лучшей позиции

    Bug(size_t dimensions, float pos_min, float pos_max, float vel_min, float vel_max) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> pos_dist(pos_min, pos_max);
        std::uniform_real_distribution<> vel_dist(vel_min, vel_max);

        for (size_t i = 0; i < dimensions; ++i) {
            position.push_back(pos_dist(gen));
            velocity.push_back(vel_dist(gen));
        }
        best_position = position;
        best_value = k::f(position); // Инициализируем значение фитнесс-функции
    }
};

// Класс роя
class Swarm {
public:
    std::vector<Bug> bugs;
    std::vector<float> global_best_position;
    float global_best_value;

    float c1, c2;     // Коэффициенты ускорения
    float inertia;    // Инерционная компонента
    size_t dimensions; // Размерность задачи
    size_t max_iterations;

    Swarm(size_t swarm_size, size_t dimensions, float pos_min, float pos_max,
          float vel_min, float vel_max, float c1, float c2, float inertia, size_t max_iterations)
        : c1(c1), c2(c2), inertia(inertia), dimensions(dimensions), max_iterations(max_iterations) {

        // Инициализируем частицы
        for (size_t i = 0; i < swarm_size; ++i) {
            bugs.emplace_back(dimensions, pos_min, pos_max, vel_min, vel_max);
        }

        // Инициализируем глобальное лучшее значение
        global_best_position = bugs[0].best_position;
        global_best_value = bugs[0].best_value;
        for (const auto& bug : bugs) {
            if (bug.best_value < global_best_value) {
                global_best_value = bug.best_value;
                global_best_position = bug.best_position;
            }
        }
    }

    Bug optimize() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.0, 1.0);

        for (size_t iter = 0; iter < max_iterations; ++iter) {
            for (auto& bug : bugs) {
                for (size_t j = 0; j < dimensions; ++j) {
                    float r1 = dist(gen);
                    float r2 = dist(gen);

                    // Корректируем скорость
                    bug.velocity[j] = inertia * bug.velocity[j]
                        + c1 * r1 * (bug.best_position[j] - bug.position[j])
                        + c2 * r2 * (global_best_position[j] - bug.position[j]);

                    // Обновляем позицию
                    bug.position[j] += bug.velocity[j];
                }

                // Вычисляем новое значение фитнесс-функции
                float fitness = k::f(bug.position);
                if (fitness < bug.best_value) {
                    bug.best_value = fitness;
                    bug.best_position = bug.position;
                }

                // Обновляем глобальное лучшее значение
                if (fitness < global_best_value) {
                    global_best_value = fitness;
                    global_best_position = bug.position;
                }
            }

            // Выводим промежуточные результаты
            if (iter % 10 == 0)
                std::cout << "Iteration " << iter  << ": Best Value = " << global_best_value << "\n";
        }

        return *std::min_element(bugs.begin(), bugs.end(), [](const Bug& a, const Bug& b)
        {
            return k::f(a.position) < k::f(b.position);
        });
    }
};

namespace k
{
inline float f(const std::vector<float> &x)
{
    float result = 10.0 * x.size();
    for (float xi : x) {
        result += xi * xi - 10.0 * std::cos(2.0 * M_PI * xi);
    }
    return result;
}

inline std::vector<float> f(const Swarm &x)
{
    std::vector<float> _;
    for (const auto &__ : x.bugs)
        _.push_back(f(__.position));
    return _;
}
}  // namespace k

inline void generateData(std::vector<std::vector<float>> &x_list, std::vector<float> &y)
{
    float step = 0.1; // Шаг между точками

    std::vector<float> x1;
    std::vector<float> x2;
    for (float i = k::MIN; i <= k::MAX; i += step) {
        x1.push_back(i);
        x2.push_back(i);
    }

    for (float xi : x1) {
        for (float xj : x2) {
            x_list[0].push_back(xi);
            x_list[1].push_back(xj);
            y.push_back(k::f({xi, xj}));
        }
    }
}

int main() {

    size_t runs = 100;

    bocchie::mark _1("1"), _2("2"), _3("3");

    size_t swarm_size = 30;
    size_t dimensions = 2;
//
//    for (auto i = 0; i < runs; ++i)
//    {
//        Swarm swarm(swarm_size,
//                    3,
//                    k::MIN,
//                    k::MAX,
//                    k::VEL_MIN,
//                    k::VEL_MAX,
//                    k::C1,
//                    k::C2,
//                    k::INERTIA,
//                    k::MAX_ITERATIONS);
//
//        auto best_bug = _1.run([&](){return swarm.optimize();});
//    }
//
//    std::cout << "For dimensions = 3:\n" << _1.to_json<bocchie::accuracy::milliseconds>() << std::endl;
//
//    for (auto i = 0; i < runs; ++i)
//    {
//        Swarm swarm(swarm_size,
//                    5,
//                    k::MIN,
//                    k::MAX,
//                    k::VEL_MIN,
//                    k::VEL_MAX,
//                    k::C1,
//                    k::C2,
//                    k::INERTIA,
//                    k::MAX_ITERATIONS);
//
//        auto best_bug = _2.run([&](){return swarm.optimize();});
//    }
//
//    std::cout << "For dimensions = 5:\n" << _2.to_json<bocchie::accuracy::milliseconds>() << std::endl;
//
//    for (auto i = 0; i < runs; ++i)
//    {
//        Swarm swarm(swarm_size,
//                    10,
//                    k::MIN,
//                    k::MAX,
//                    k::VEL_MIN,
//                    k::VEL_MAX,
//                    k::C1,
//                    k::C2,
//                    k::INERTIA,
//                    k::MAX_ITERATIONS);
//
//        auto best_bug = _3.run([&](){return swarm.optimize();});
//    }
//
//    std::cout << "For dimensions = 10:\n" << _3.to_json<bocchie::accuracy::milliseconds>() << std::endl;


    Swarm swarm(swarm_size,
                dimensions,
                k::MIN,
                k::MAX,
                k::VEL_MIN,
                k::VEL_MAX,
                k::C1,
                k::C2,
                k::INERTIA,
                k::MAX_ITERATIONS);


    auto print_swarm_plot = [&](){
        std::vector<std::vector<float>> x_list(dimensions);
        std::vector<float> y;

        generateData(x_list, y);
        std::vector<std::vector<float>> p_list;
        p_list.reserve(swarm_size);
        for (const auto &_p : swarm.bugs) {
            p_list.emplace_back();
            for (const auto &gene : _p.position)
                p_list.back().push_back(gene);
        }

        std::vector<float> p;
        for (const auto &_p : swarm.bugs)
            p.push_back(k::f(_p.position));

        plot::print3d(x_list, y, p_list, p);
    };

    print_swarm_plot();

    auto best_bug = swarm.optimize();

    print_swarm_plot();


    {
        std::vector<std::vector<float>> x_list(dimensions);
        std::vector<float> y;

        generateData(x_list, y);
        std::vector<std::vector<float>> p_list;
        p_list.push_back(best_bug.position);


        std::vector<float> p;
        p.push_back(k::f(best_bug.position));

        plot::print3d(x_list, y, p_list, p);
    }

    return 0;
}

}

#endif //LAB7_HPP
