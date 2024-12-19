#ifndef LAB1
#define LAB1

#include <iostream>
#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <exception>
#include <random>
#include "Plot.h"

namespace lab1
{

typedef uint16_t Entity;
typedef std::vector<Entity> Population;
typedef std::vector<float> Fitness;

namespace k
{
const uint32_t PRECISION = 1 << 14;

const uint8_t START = 1, END = 10;

float P_MUTATION = 0.0001;  // Вероятность мутации
float P_CROSSING = 0.5;     // Вероятность скрещивания

uint16_t N = 100;

float f(float x)
{
    return std::log(x) * std::cos(3 * x - 15);
}

template<typename T>
std::vector<float> f(std::vector<T> _x)
{
    std::vector<float> __x;
    __x.reserve(_x.size());
    for (const auto &x : _x)
        __x.push_back(f(float(x)));
    return __x;
}

float fitness(float x)
{
    return - std::log(x) * std::cos(3 * x - 15);
}

template<typename T>
Fitness fitness(std::vector<T> _x)
{
    std::vector<float> __x;
    __x.reserve(_x.size());
    for (const auto &x : _x)
        __x.push_back(fitness(float(x)));
    return __x;
}
// Функция для перевода вещественного числа из промежутка в целочисленный эквивалент
uint16_t chunk(float value)
{
    float step = (float)(k::END - k::START) / k::PRECISION;
    return uint16_t((value - k::START) / float(step));
}

// Функция для обратного оперевода
float value(uint16_t chunk)
{
    float step = (float)(k::END - k::START) / k::PRECISION;
    return k::START + (float)chunk * step;
}
}

template<typename T1, typename T2>
void is_aboba(std::vector<T1> _, std::vector<T2> __)
{
    if (_.size() != __.size())
        throw std::exception();
}

namespace op
{
Population selection(const Population &population,
                     const Fitness &fitness,
                     int num_selected)
{
    is_aboba(population, fitness);

    std::vector<std::pair<float /* <- fit */ , uint16_t /* <- source ind */>> sorted_fitness;
    for (size_t i = 0; i < fitness.size(); ++i)
        sorted_fitness.emplace_back(fitness[i], i);

//    auto min = *std::min_element(fitness.begin(), fitness.end());
//
//    auto positive_fitness = fitness;
//    for (auto &fit : positive_fitness)
//        fit += std::abs(min);
//
//    auto sum = std::accumulate(positive_fitness.begin(), positive_fitness.end(), 0.0f);
//
//    std::vector<std::pair<float /* <- prob */ , uint16_t /* <- source ind */>> probabilities;
//    for (size_t i = 0; i < positive_fitness.size(); ++i)
//        probabilities.push_back({positive_fitness[i] / float(sum), i});

    std::sort(sorted_fitness.begin(), sorted_fitness.end(), [](const auto& a, const auto& b)
    {
        return a.first > b.first;
    });

    Population selected_population;
    for (int i = 0; i < num_selected; ++i) {
        selected_population.push_back(population[sorted_fitness[i].second]);
    }

    return selected_population;
}

Population crossover(const Population &population)
{
    Population new_population;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(1, 15);

    for (size_t i = 0; i < population.size(); i += 2) {
        uint16_t parent1 = population[i];
        uint16_t parent2 = (i + 1 < population.size()) ? population[i + 1] : population[i];

        int crossover_point = dist(gen);

        uint16_t mask1 = (1 << crossover_point) - 1;
        uint16_t mask2 = ~mask1;

        uint16_t child1 = (parent1 & mask1) | (parent2 & mask2);

        new_population.push_back(child1);
    }

    return new_population;
}

Population mutate(const Population &population, float mutation_rate)
{
    Population mutated_population;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (uint16_t entity : population) {
        for (int i = 0; i < 14; ++i) {
            if (dist(gen) < mutation_rate) {
                entity ^= (1 << i);  // Инвертируем i-й бит
            }
        }

        // Ограничиваем значение в пределах 2**14
        if (entity > (1 << 14) - 1) {
            entity = (1 << 14) - 1;
        }

        mutated_population.push_back(entity);
    }

    return mutated_population;
}
}
int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(1, 10);

    Population population;
    population.reserve(k::N);
    for (size_t i = 0; i < k::N; ++i) {
        population.push_back(k::chunk((float)dist(gen)));
    }

    Population source_population = population;


    std::vector<float> x(k::N);
    double step = float(k::END - k::START) / k::N;
    double current_value = 1.0;
    for (int i = 0; i < k::N; ++i) {
        x[i] = current_value;
        current_value += step;
    }

    std::vector<float> y;
    for (auto _ : x)
        y.push_back(k::f(_));

    {
        std::vector<float> _x;
        for (const auto& _p: population)
            _x.push_back(k::value(_p));

        std::vector<float> p;
        for (const auto &_p : population)
            p.push_back(k::f(k::value(_p)));


        plot::print(x, y, _x, p);
    }

    auto get_unchunked_population = [](const Population& _){
        std::vector<float> out;
        out.reserve(_.size());
        for (const auto& chunk: _)
            out.push_back(k::value(chunk));
        return out;
    };

    for (int i = 0; i < 50; ++i)
    {
        // 1. Выбор родителей для процесса размножения (работает оператор селекции - репродукции)
        Population parents = op::selection(population, k::fitness(get_unchunked_population(population)), k::N * k::P_CROSSING);

        // 2. Создание потомков выбранных пар родителей (работает оператор скрещивания - кроссинговера)
        Population childs = op::crossover(parents);

        // 3. Мутация новых особей (работает оператор мутации)
        childs = op::mutate(childs, k::P_MUTATION);

        // 4. Расширение популяции за счет добавления новых только что порожденных особей
        population.insert(population.end(), childs.begin(), childs.end());

        // 5. Сокращение расширенной популяции до исходного размера (работает оператор редукции)
        population = op::selection(population, k::fitness(get_unchunked_population(population)), k::N);
    }

    {
        std::vector<float> _x;
        for (const auto& _p: population)
            _x.push_back(k::value(_p));

        std::vector<float> p;
        for (const auto &_p : population)
            p.push_back(k::f(k::value(_p)));

        plot::print(x, y, _x, p);
    }

    return 0;
}


}

#endif