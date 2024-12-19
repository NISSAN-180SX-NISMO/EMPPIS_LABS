#ifndef LAB5_HPP
#define LAB5_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <cstdint>
#include "Plot.h"
#include "Random.hpp"
#include "Bocchiemark.hpp"

namespace lab5
{

constexpr uint8_t n = 2;


struct Entity
{
    std::array<float, n> values;
    std::array<float, n> biases;
};

typedef std::vector<Entity> Population;

typedef std::vector<float> Fitness;

namespace k
{
const size_t EPOCH_CRITERIA = 300;
const size_t EPOCH_LIMIT = 200'000;
const size_t FLOAT_EQ_PRECISION = 10;

const float MIN = -5.12, MAX = 5.12;

struct P
{
    float _CROSSING = 0.5;     // Вероятность скрещивания
    float _MUTATION = 0.0001;  // Вероятность мутации
};

uint16_t N = 100;

inline float f(const Entity &x)
{
    float result = 10.0 * x.values.size();
    for (float xi : x.values) {
        result += xi * xi - 10.0 * std::cos(2.0 * M_PI * xi);
    }
    return result;
}

inline std::vector<float> f(const Population &x)
{
    std::vector<float> _;
    for (const auto &__ : x)
        _.push_back(f(__));
    return _;
}

// Функция для сравнения двух чисел с плавающей запятой с точностью до n знаков
bool are_equal_with_precision(float a, float b, int precision = FLOAT_EQ_PRECISION) {
    // Степень 10 для проверки точности
    float epsilon = std::pow(10.0f, -precision);

    // Проверяем разницу между числами
    return std::fabs(a - b) < epsilon;
}

void applyReshenbergRule(Population& population, size_t success, size_t total) {
    double successRatio = static_cast<double>(success) / total;

    // Если доля успеха больше 1/5
    if (successRatio > 1.0 / 5.0) {
        for (auto& individual : population) {
            for (auto& param : individual.biases) {
                param /= 0.82;  // уменьшение на 18%
            }
        }
    }
        // Если доля успеха меньше 1/5
    else if (successRatio < 1.0 / 5.0) {
        for (auto& individual : population) {
            for (auto& param : individual.biases) {
                param *= 0.82;  // увеличение на 18%
            }
        }
    }
}

}

namespace op
{
Population selection(const Population &population,
                     const Fitness &fitness,
                     int num_selected)
{
    std::vector<std::pair<float /* <- fit */ , uint16_t /* <- source ind */>> sorted_fitness;
    for (size_t i = 0; i < fitness.size(); ++i)
        sorted_fitness.emplace_back(fitness[i], i);


    std::sort(sorted_fitness.begin(), sorted_fitness.end(), [](const auto &a, const auto &b)
    {
        return a.first < b.first;
    });

    Population selected_population;
    for (int i = 0; i < num_selected; ++i) {
        selected_population.push_back(population[sorted_fitness[i].second]);
    }

    return selected_population;
}

Population crossover(const Population &parents)
{
    Population children;

    // Размер родителей
    size_t n = parents[0].values.size();

    // Для каждой пары родителей создаем одного потомка
    for (size_t i = 0; i < parents.size(); i += 2) {
        Entity child1, child2;
        const Entity &parent1 = parents[i];
        const Entity &parent2 = parents[i + 1];

        // Применяем однородный кроссовер для каждого гена
        for (size_t j = 0; j < n; ++j) {
            if (random::_float(0.f, 1.f) < 0.5) {
                child1.values[j] = parent1.values[j];
                child2.values[j] = parent2.values[j];
            } else {
                child1.values[j] = parent2.values[j];
                child2.values[j] = parent1.values[j];
            }
        }

        children.push_back(child1);
        children.push_back(child2);
    }

    return children;
}

Population mutate(const Population &population, float mutation_rate, size_t &successful_mutations)
{
    Population mutated_population = population;

    for (size_t i = 0; i < population.size(); ++i) {
        Entity &entity = mutated_population[i];
        const Entity &original_entity = population[i];

        if (random::_float(0.f, 1.f) < mutation_rate) {
            // Сохраняем исходное значение для отката, если мутация не успешна
            decltype(n) lucker = random::_int(0, n - 1);
            float original_value = entity.values[lucker];

            // Генерируем мутацию
            float mutation = random::_float(k::MIN, k::MAX);
            entity.values[lucker] = mutation;

            // Проверяем успешность мутации
            if (k::f(entity) < k::f(original_entity)) {
                ++successful_mutations;
            } else {
                // Откатываем мутацию, если она не успешна
                entity.values[lucker] = original_value;
            }
        }
    }

    return mutated_population;
}

}

// Генерация данных для построения графика
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

struct aboba
{
    Entity x;
    float y;
    std::string to_string()
    {
        std::string out;
        for (size_t i = 0; i < n; ++i)
            out += "x" + std::to_string(i) + " = " + std::to_string(x.values[i]) + ", ";

        out += "y = " + std::to_string(y);
        return out;
    }
};

Population genRandPopulation()
{
    Population population;
    population.reserve(k::N);
    for (size_t i = 0; i < k::N; ++i) {
        Entity _;
        for (size_t j = 0; j < n; ++j)
        {
            _.values[j] = random::_float(k::MIN, k::MAX);
            _.biases[j] = random::_float(0, 1);
        }
        population.push_back(_);
    }
    return population;
}

std::pair<int, aboba> _main(k::P p, Population population, bool show_plots, std::string bmark_name = "")
{
    Population source_population = population;

    std::vector<std::vector<float>> x_list(n);
    std::vector<float> y;

    generateData(x_list, y);

    if (show_plots) {
        std::vector<std::vector<float>> p_list;
        p_list.reserve(population.size());
        for (const auto &_p : population) {
            p_list.emplace_back();
            for (const auto &gene : _p.values)
                p_list.back().push_back(gene);
        }

        std::vector<float> p;
        for (const auto &_p : population)
            p.push_back(k::f(_p));

        plot::print3d(x_list, y, p_list, p);
    }

    bool exit = false;
    size_t epoch_counter = 0;

    auto min = [](const Population& population)
    {
        return *std::min_element(population.begin(), population.end(),
                                 [](const auto &a, const auto &b)
                                 {
                                     return k::f(a) < k::f(b);
                                 });
    };

    auto fine_optimum = [&epoch_counter, &min](const Population &population) -> bool
    {
        static size_t best_optimum_counter = 0;
        static float last_optimum = 0xFFFFFF;

        float optimum = k::f(min(population));
        if (not k::are_equal_with_precision(last_optimum, optimum))
        {
            last_optimum = optimum;
            best_optimum_counter = 0;
        }
        else best_optimum_counter++;

        return best_optimum_counter == k::EPOCH_CRITERIA || epoch_counter == k::EPOCH_LIMIT;
    };


//    for (int i = 0; i < 50; ++i)

    size_t successful_mutations = 0;

    while (not exit) {
        epoch_counter++;

        // 1. Выбор родителей для процесса размножения (работает оператор селекции - репродукции)
        Population parents = op::selection(population, k::f(population), k::N * p._CROSSING);

        // 2. Создание потомков выбранных пар родителей (работает оператор скрещивания - кроссинговера)
        Population childs = op::crossover(parents);

        // 3. Мутация новых особей (работает оператор мутации)
        childs = op::mutate(childs, p._MUTATION, successful_mutations);

        // 4. Расширение популяции за счет добавления новых только что порожденных особей
        population.insert(population.end(), childs.begin(), childs.end());

        // 5. Сокращение расширенной популяции до исходного размера (работает оператор редукции)
        population = op::selection(population, k::f(population), k::N);

        k::applyReshenbergRule(population, successful_mutations, epoch_counter);

        exit = fine_optimum(population);
    }


    if (show_plots) {
        std::vector<std::vector<float>> p_list;
        p_list.reserve(population.size());
        for (const auto &_p : population) {
            p_list.emplace_back();
            for (const auto &gene : _p.values)
                p_list.back().push_back(gene);
        }

        std::vector<float> p;
        for (const auto &_p : population)
            p.push_back(k::f(_p));

        plot::print3d(x_list, y, p_list, p);
    }

    return {(int) epoch_counter, aboba { min(population), k::f(min(population)) } };
}


namespace stats {

// Функция для нахождения минимального значения
int findMin(const std::vector<int>& vec) {
    return *std::min_element(vec.begin(), vec.end());
}

aboba findMin(const std::vector<aboba>& vec) {
    return *std::min_element(vec.begin(), vec.end(), [](const auto& a, const auto& b)
    {
        return a.y < b.y;
    });
}

// Функция для нахождения максимального значения
int findMax(const std::vector<int>& vec) {
    return *std::max_element(vec.begin(), vec.end());
}

// Функция для нахождения среднего значения
double findMean(const std::vector<int>& vec) {
    double sum = std::accumulate(vec.begin(), vec.end(), 0);
    return sum / vec.size();
}

double findMean(const std::vector<float>& vec) {
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0f);
    return sum / vec.size();
}

// Функция для нахождения медианы
double findMedian(const std::vector<int>& vec) {
    auto copy = vec;
    size_t size = copy.size();
    std::sort(copy.begin(), copy.end());

    if (size % 2 == 0) {
        return (copy[size / 2 - 1] + copy[size / 2]) / 2.0;
    } else {
        return copy[size / 2];
    }
}

aboba findMedian(const std::vector<aboba>& vec) {
    auto copy = vec;
    size_t size = copy.size();
    std::sort(copy.begin(), copy.end(), [](const auto& a, const auto& b)
    {
        return a.y < b.y;
    });

    return copy[size / 2];
}
}

inline int main()
{
    size_t runs = 1;

    std::vector<Population> populations;
    for (size_t i = 0; i < runs; ++i)
        populations.push_back(genRandPopulation());

    constexpr bocchie::accuracy mark_accuracy = bocchie::accuracy::milliseconds;
    bocchie::mark _1("1"), _2("2"), _3("3"), _4("4");

    auto bocchiemark = [&](float crossing_prob, float mutation_prob, const bocchie::mark _mark){
        std::shared_ptr<bocchie::mark> mark = std::make_shared<bocchie::mark>(_mark);
        std::cout << "-----------------------------------------------\n"
                  << "Start education with Crossing prob = " << crossing_prob << " and Mutation prob = " << mutation_prob << " ...\n";

        std::vector<int> epochs; epochs.reserve(runs);
        std::vector<aboba> optimums; optimums.reserve(runs);
        for (size_t i = 0; i < runs; ++i)
        {
            auto [epoch, aboba] = mark->run(_main, k::P{crossing_prob, mutation_prob}, populations[i], true, mark->get_runnable_name());
            epochs.push_back(epoch);
            optimums.push_back(aboba);
        }

        std::cout << "Min epochs: " << stats::findMin(epochs) << std::endl;
        std::cout << "Max epochs: " << stats::findMax(epochs) << std::endl;
        std::cout << "Avg epochs: " << stats::findMean(epochs) << std::endl;
        std::cout << "Median epochs: " << stats::findMedian(epochs) << std::endl;

        std::cout << "MIN OPTIMUM: " << stats::findMin(optimums).to_string() << std::endl;
        std::cout << "MEDIAN OPTIMUM: " << stats::findMedian(optimums).to_string() << std::endl;

        std::cout << mark->to_json<mark_accuracy>() << std::endl;
    };

//    bocchiemark(0.5, 0.0001, _1);
//    bocchiemark(0.2, 0.0001, _2);
//    bocchiemark(0.7, 0.0001, _3);
    bocchiemark(0.2, 0.1, _4);

    return 0;
}
}

#endif //LAB5_HPP
