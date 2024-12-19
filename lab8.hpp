#ifndef LAB8_HPP
#define LAB8_HPP

#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "Random.hpp"

namespace lab8
{
bocchie::mark selection_mark("selection"), crossover_mark("crossover"), mutate_mark("mutate"),
    insert_mark("population.insert"), EPOCH_MARK("EPOCH"), list_fitness_mark("list fitness");

struct Entity
{
    double a, b;
    [[nodiscard]] double evaluate_Ef(double L) const
    {
        return a + std::pow(L, b);
    }
} bestEntity;

typedef std::vector<Entity> Population;

struct Project
{
    struct Feature
    {
        double L;  // В килостроках
    } feature;
    struct Target
    {
        double Ef; // Реальная стоимость
    } target;
};

namespace dataset
{

std::vector<Project> projects_dataset;

std::vector<Project> train_set;

std::vector<Project> test_set;

}  // namespace dataset

namespace k
{

constexpr uint16_t N = 1000;

struct P
{
    double c, m;
};

constexpr double MUTATION_RANGE = 0.1;

constexpr double MUTATION_OFFSET = 0.05;

namespace cost
{

double MMRE(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size() || a.empty()) {
        throw std::invalid_argument("Vectors sizes not equal or empty =(");
    }

    double totalRelativeError = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        totalRelativeError += std::abs((a[i] - b[i]) / a[i]);
    }

    return totalRelativeError / (double) a.size();
}

}  // namespace cost

double fitness(const Entity &entity, const std::vector<Project> &set = dataset::train_set)
{
    std::vector<double> predictions, reality;
    predictions.reserve(set.size());
    reality.reserve(set.size());

    for (const auto &project : set) {
        predictions.push_back(entity.evaluate_Ef(project.feature.L));
        reality.push_back(project.target.Ef);
    }

    return cost::MMRE(predictions, reality);
}

std::vector<double> fitness(const Population &population, const std::vector<Project> &set = dataset::train_set)
{
    std::vector<double> out;
    out.reserve(population.size());

    for (const auto &entity : population)
        out.emplace_back(fitness(entity, set));

    return out;
}

}  // namespace k

namespace op
{

// Ранговый оператор селекции
Population selection(const Population &population, const std::vector<double> &fitness, double crossing_rate)
{
    Population selected_population;
    std::vector<std::pair<double, Entity>> ranked_population;

    // Ранжируем популяцию по фитнес-функции
    for (size_t i = 0; i < population.size(); ++i) {
        ranked_population.emplace_back(fitness[i], population[i]);
    }

    // Сортируем популяцию по фитнес-функции
    std::sort(ranked_population.begin(), ranked_population.end(), [](const auto &a, const auto &b)
    {
        return a.first < b.first;
    });

    // Выбираем лучшие особи
    for (size_t i = 0; i < size_t(population.size() * crossing_rate); ++i) {
        selected_population.push_back(ranked_population[i].second);
    }

    bestEntity = selected_population.front();
    return selected_population;
}

// Арифметический оператор скрещивания
Population crossover(const Population &population)
{
    Population new_population;

    for (size_t i = 0; i < population.size(); i += 2) {
        if (i + 1 < population.size()) {
            const Entity &parent1 = population[i];
            const Entity &parent2 = population[i + 1];

            Entity child1, child2;
            double alpha = random::_double(0.0, 1.0);

            child1.a = alpha * parent1.a + (1 - alpha) * parent2.a;
            child1.b = alpha * parent1.b + (1 - alpha) * parent2.b;

            child2.a = (1 - alpha) * parent1.a + alpha * parent2.a;
            child2.b = (1 - alpha) * parent1.b + alpha * parent2.b;

            new_population.push_back(child1);
            new_population.push_back(child2);
        }
        else {
            new_population.push_back(population[i]);
        }
    }

    return new_population;
}

// Арифметический оператор мутации
Population mutate(const Population &population, double mutation_rate)
{
    Population mutated_population = population;

    for (auto &entity : mutated_population) {
        if (random::_double(0.0, 1.0) < mutation_rate) {
            entity.a +=
                random::_double(0.0, 1.0) * k::MUTATION_RANGE - k::MUTATION_OFFSET; // Мутация с небольшим изменением
            entity.b +=
                random::_double(0.0, 1.0) * k::MUTATION_RANGE - k::MUTATION_OFFSET; // Мутация с небольшим изменением
        }
    }

    return mutated_population;
}

}  // namespace op

Entity _main(Population population, k::P p, size_t epochs = 1000, size_t log_epoch = 100)
{
    for (size_t epoch = 0; epoch < epochs; ++epoch)
    {
        EPOCH_MARK.run([&]()
       {
           Population parents = selection_mark.run(op::selection, population, list_fitness_mark.run([&](){return k::fitness(population);}), p.c);

           Population childs = crossover_mark.run(op::crossover,parents);

           childs = mutate_mark.run(op::mutate,childs, p.m);

           insert_mark.run([&](){return population.insert(population.end(), childs.begin(), childs.end());});

           population = selection_mark.run(op::selection, population, list_fitness_mark.run([&](){return k::fitness(population);}), (double)k::N / population.size());

           if (log_epoch != 0 && epoch % log_epoch == 0)
               std::cout << "Epoch " << epoch << ": best fitness (MMRE) = " << list_fitness_mark.run([&](){return k::fitness(bestEntity);}) << std::endl;

           return true;
       });
    }

    return bestEntity;
}

void init_dataset()
{
    // Данные из таблицы 8.2
    {
        using namespace dataset;
        projects_dataset.emplace_back(Project::Feature{2.2}, Project::Target{8.4});
        projects_dataset.emplace_back(Project::Feature{3.5}, Project::Target{10.8});
        projects_dataset.emplace_back(Project::Feature{5.5}, Project::Target{18});
        projects_dataset.emplace_back(Project::Feature{6.0}, Project::Target{24});
        projects_dataset.emplace_back(Project::Feature{9.7}, Project::Target{25.2});
        projects_dataset.emplace_back(Project::Feature{7.7}, Project::Target{31.2});
        projects_dataset.emplace_back(Project::Feature{11.3}, Project::Target{36});
        projects_dataset.emplace_back(Project::Feature{8.2}, Project::Target{36});
        projects_dataset.emplace_back(Project::Feature{6.5}, Project::Target{42});
        projects_dataset.emplace_back(Project::Feature{8.0}, Project::Target{42});
        projects_dataset.emplace_back(Project::Feature{20.0}, Project::Target{48});
        projects_dataset.emplace_back(Project::Feature{10.0}, Project::Target{48});
        projects_dataset.emplace_back(Project::Feature{15.0}, Project::Target{48});
        projects_dataset.emplace_back(Project::Feature{10.4}, Project::Target{50});
        projects_dataset.emplace_back(Project::Feature{13.0}, Project::Target{60});
        projects_dataset.emplace_back(Project::Feature{14.0}, Project::Target{60});
        projects_dataset.emplace_back(Project::Feature{19.7}, Project::Target{60});
        projects_dataset.emplace_back(Project::Feature{32.5}, Project::Target{60});
        projects_dataset.emplace_back(Project::Feature{31.5}, Project::Target{60});
        projects_dataset.emplace_back(Project::Feature{12.5}, Project::Target{62});
        projects_dataset.emplace_back(Project::Feature{15.4}, Project::Target{70});
        projects_dataset.emplace_back(Project::Feature{20.0}, Project::Target{72});
        projects_dataset.emplace_back(Project::Feature{7.5}, Project::Target{72});
        projects_dataset.emplace_back(Project::Feature{16.3}, Project::Target{82});
        projects_dataset.emplace_back(Project::Feature{15.0}, Project::Target{90});
        projects_dataset.emplace_back(Project::Feature{11.4}, Project::Target{98.8});
        projects_dataset.emplace_back(Project::Feature{21.0}, Project::Target{107});
        projects_dataset.emplace_back(Project::Feature{16.0}, Project::Target{114});
        projects_dataset.emplace_back(Project::Feature{25.9}, Project::Target{117.6});
        projects_dataset.emplace_back(Project::Feature{24.6}, Project::Target{117.6});
        projects_dataset.emplace_back(Project::Feature{29.5}, Project::Target{120});
        projects_dataset.emplace_back(Project::Feature{19.3}, Project::Target{155});
        projects_dataset.emplace_back(Project::Feature{32.6}, Project::Target{170});
        projects_dataset.emplace_back(Project::Feature{35.5}, Project::Target{192});
        projects_dataset.emplace_back(Project::Feature{38.0}, Project::Target{210});
        projects_dataset.emplace_back(Project::Feature{48.5}, Project::Target{239});
        projects_dataset.emplace_back(Project::Feature{47.5}, Project::Target{252});
        projects_dataset.emplace_back(Project::Feature{70.0}, Project::Target{278});
        projects_dataset.emplace_back(Project::Feature{66.6}, Project::Target{300});
        projects_dataset.emplace_back(Project::Feature{66.6}, Project::Target{352.8});
        projects_dataset.emplace_back(Project::Feature{50.0}, Project::Target{370});
        projects_dataset.emplace_back(Project::Feature{79.0}, Project::Target{400});
        projects_dataset.emplace_back(Project::Feature{90.0}, Project::Target{450});
        projects_dataset.emplace_back(Project::Feature{78.0}, Project::Target{571.4});
        projects_dataset.emplace_back(Project::Feature{100.0}, Project::Target{215});
        projects_dataset.emplace_back(Project::Feature{150.0}, Project::Target{324});
        projects_dataset.emplace_back(Project::Feature{100.0}, Project::Target{360});
        projects_dataset.emplace_back(Project::Feature{100.0}, Project::Target{360});
        projects_dataset.emplace_back(Project::Feature{190.0}, Project::Target{420});
        projects_dataset.emplace_back(Project::Feature{115.8}, Project::Target{480});

        // Перемешиваем данные для случайного разделения
        std::srand(static_cast<unsigned>(std::time(nullptr)));
        std::shuffle(projects_dataset.begin(), projects_dataset.end(), std::mt19937(std::random_device()()));

        // Разделяем на обучающую и тестовую выборки
        size_t test_size = projects_dataset.size() - 40;
        size_t train_size = projects_dataset.size() - test_size;

        train_set.insert(train_set.end(), projects_dataset.begin(), projects_dataset.begin() + train_size);
        test_set.insert(test_set.end(), projects_dataset.begin() + train_size, projects_dataset.end());
    };
}

void main()
{
    init_dataset();

    Population population = []()
    {
        Population population;
        population.reserve(k::N);
        for (auto i = 0; i < k::N; ++i)
            population.emplace_back(random::_double(0, 2), random::_double(0, 2));
        return population;
    }();

    auto best = _main(population, k::P{0.7, 0.1}, 3000, 100);

    std::cout << selection_mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;
    std::cout << crossover_mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;
    std::cout << mutate_mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;
    std::cout << insert_mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;
    std::cout << list_fitness_mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;
    std::cout << EPOCH_MARK.to_json<bocchie::accuracy::milliseconds>() << std::endl;

    std::cout << "Best entity have a = " << best.a << ", b = " << best.b << ".\n"
              << "Check on test set:\n";

    for (const auto &project : dataset::test_set)
        std::cout << "L = "<< std::setw(7) << project.feature.L
                  << "\tPredicted value = "<< std::setw(7) << best.evaluate_Ef(project.feature.L)
                  << "\tReality value = " << std::setw(7) << project.target.Ef << std::endl;

    std::cout << "Fitness (MMRE) on test set = " << k::fitness(best, dataset::test_set) << std::endl;
}

}  // namespace lab8

#endif //LAB8_HPP
