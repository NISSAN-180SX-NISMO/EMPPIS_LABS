#ifndef LAB3_HPP
#define LAB3_HPP
#include "Random.hpp"

namespace lab3
{
typedef std::pair<float, float> City;
typedef std::vector<City> Tour;
typedef std::vector<float> Fitness;

struct Entity
{
    Tour tour;
};

typedef std::vector<Entity> Population;

namespace k
{
const size_t EPOCH_CRITERIA = 1000;
const size_t EPOCH_LIMIT = 200'000;
const size_t FLOAT_EQ_PRECISION = 10;

const uint16_t N = 1000;

struct P
{
    float _CROSSING = 0.5;     // Вероятность скрещивания
    float _MUTATION = 0.0001;  // Вероятность мутации
};


// Функция для вычисления евклидова расстояния между двумя точками
float calculateDistance(const City& a, const City& b) {
    float dx = a.first - b.first;
    float dy = a.second - b.second;
    return std::sqrt(dx * dx + dy * dy);
}

// Вычисление фитнеса для одного Entity (общая длина маршрута)
float fitness(const Tour& tour) {
    if (tour.size() < 2) return 0.0f; // Если меньше двух городов, длина маршрута равна 0

    float totalDistance = 0.0f;

    for (size_t i = 0; i < tour.size(); ++i) {
        const City& from = tour[i];
        const City& to = tour[(i + 1) % tour.size()]; // Кольцевой маршрут
        totalDistance += calculateDistance(from, to);
    }

    return totalDistance;
}

// Вычисление фитнеса для всей популяции
std::vector<float> fitness(const Population& population) {
    std::vector<float> fitnessValues;
    fitnessValues.reserve(population.size()); // Резервируем место для оптимизации

    for (const Entity& entity : population) {
        fitnessValues.push_back(fitness(entity.tour)); // Вычисляем фитнес каждого Entity
    }

    return fitnessValues;
}

// Функция для сравнения двух чисел с плавающей запятой с точностью до n знаков
bool are_equal_with_precision(float a, float b, int precision = FLOAT_EQ_PRECISION) {
    return a == b;
    // Степень 10 для проверки точности
    float epsilon = std::pow(10.0f, -precision);

    // Проверяем разницу между числами
    return std::fabs(a - b) < epsilon;
}

}  // namespace k

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

Population crossover(const Population &population)
{
    Population offspring;

    for (size_t i = 0; i < population.size(); i++)
    {
        // Выбираем двух случайных родителей
        const auto &parent1 = population[random::_int(0, population.size() - 1)].tour;
        const auto &parent2 = population[random::_int(0, population.size() - 1)].tour;

        // Создаём потомка, равного размеру родителей
        Tour child(parent1.size());
        std::fill(child.begin(), child.end(), City{-1, -1});
        child.front() = parent1.front();
        child.back() = parent1.back();

        // Выбираем случайный промежуток для обмена, исключая первый и последний город
        size_t start = random::_int(1, parent1.size() - 2); // Исключаем первый и последний города
        size_t end = random::_int(start, parent1.size() - 2); // Исключаем первый и последний города

        if (start > end)
            std::swap(start, end);

        // Копируем часть от первого родителя
        for (size_t j = start; j <= end; ++j)
        {
            child[j] = parent1[j];
        }

        // Заполняем оставшиеся города из второго родителя
        for (size_t j = 1; j < parent2.size() - 1; ++j) // Исключаем первый и последний города
        {
            if (std::find(child.begin(), child.end(), parent2[j]) == child.end())
            {
                auto it = std::find(child.begin(), child.end(), City{-1, -1});
                if (it != child.end())
                {
                    *it = parent2[j];
                }
            }
        }

        offspring.push_back({child});
    }

    return offspring;
}


// Мутация методом перестановки (swap mutation)
Population mutate(const Population &population, float mutation_rate)
{
    Population mutated_population = population;


    for (auto &entity : mutated_population)
    {
        if (random::_float(0, 1) < mutation_rate)
        {
            size_t idx1 = random::_int(1, entity.tour.size() - 2);
            size_t idx2 = random::_int(1, entity.tour.size() - 2);

            std::swap(entity.tour[idx1], entity.tour[idx2]);
        }
    }

    return mutated_population;
}

}  // namespace op

Entity createRandomEntity(std::vector<City> cityPool)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());

    // Перемешиваем города
    std::shuffle(cityPool.begin(), cityPool.end(), gen);

    cityPool.push_back(cityPool.front());

    return Entity { cityPool };
}

Entity createRandomEntityWithFixedStart(const std::vector<City>& cityPool)
{
    if (cityPool.empty())
        throw std::invalid_argument("City pool cannot be empty.");

    static std::random_device rd;
    static std::mt19937 gen(rd());

    // Извлекаем первый город
    City fixedCity = cityPool.front();

    // Копируем города (кроме первого) для случайного перемешивания
    std::vector<City> shuffledCities(cityPool.begin() + 1, cityPool.end());

    // Перемешиваем оставшиеся города
    std::shuffle(shuffledCities.begin(), shuffledCities.end(), gen);

    // Формируем тур
    std::vector<City> tour = {fixedCity};
    tour.insert(tour.end(), shuffledCities.begin(), shuffledCities.end());
    tour.push_back(fixedCity);

    return Entity{tour};
}

Population createRandomPopulation(const std::vector<City>& cityPool, bool fixedStart = true)
{
    Population out;
    for (size_t i = 0; i < k::N; ++i)
        out.push_back(fixedStart ? createRandomEntityWithFixedStart(cityPool) : createRandomEntity(cityPool));
    return out;
};

void printTour(const Tour& cityPool, const Tour& tour) {
    for (const auto& city : tour) {
        // Ищем индекс города в пуле
        auto it = std::find(cityPool.begin(), cityPool.end(), city);
        if (it != cityPool.end()) {
            // Индекс найден
            int index = std::distance(cityPool.begin(), it);
            std::cout << index + 1 << "-";
        }
    }
    std::cout << std::endl;  // Печатаем новую строку после завершения вывода
}

auto min = [](const Population& population)
{
    return *std::min_element(population.begin(), population.end(),
                             [](const auto &a, const auto &b)
                             {
                                 return k::fitness(a.tour) < k::fitness(b.tour);
                             });
};

Entity _main(Population population, k::P p)
{
    bool exit = false;
    size_t epoch_counter = 0;

    auto fine_optimum = [&epoch_counter](const Population &population) -> bool
    {
        static size_t best_optimum_counter = 0;
        static float last_optimum = 0xFFFFFF;

        float optimum = k::fitness(min(population).tour);
        if (not k::are_equal_with_precision(last_optimum, optimum))
        {
            last_optimum = optimum;
            best_optimum_counter = 0;
        }
        else best_optimum_counter++;

        return best_optimum_counter == k::EPOCH_CRITERIA || epoch_counter == k::EPOCH_LIMIT;
    };

    for (size_t i = 0; i < 2'000; ++i)
//    while (not exit)
    {
        epoch_counter++;

        // 1. Выбор родителей для процесса размножения (работает оператор селекции - репродукции)
        Population parents = op::selection(population, k::fitness(population), k::N * p._CROSSING);

        // 2. Создание потомков выбранных пар родителей (работает оператор скрещивания - кроссинговера)
        Population childs = op::crossover(parents);

        // 3. Мутация новых особей (работает оператор мутации)
        childs = op::mutate(childs, p._MUTATION);

        // 4. Расширение популяции за счет добавления новых только что порожденных особей
        population.insert(population.end(), childs.begin(), childs.end());

        // 5. Сокращение расширенной популяции до исходного размера (работает оператор редукции)
        population = op::selection(population, k::fitness(population), k::N);

        if (i % 1000 == 0)
        {
            std::cout << "Epoch " << i << ", tour length: " << k::fitness(min(population).tour) << std::endl;
//            printTour(berlin52, min(population).tour);
            plot::city_graph(min(population).tour);
        }

//        exit = fine_optimum(population);
    }

    std::cout << "Educated by " << epoch_counter << " epochs\n";

    return min(population);
}

void main()
{
    std::vector<City> berlin52;
    {
        berlin52.emplace_back(565.0, 575.0);
        berlin52.emplace_back(25.0, 185.0);
        berlin52.emplace_back(345.0, 750.0);
        berlin52.emplace_back(945.0, 685.0);
        berlin52.emplace_back(845.0, 655.0);
        berlin52.emplace_back(880.0, 660.0);
        berlin52.emplace_back(25.0, 230.0);
        berlin52.emplace_back(525.0, 1000.0);
        berlin52.emplace_back(580.0, 1175.0);
        berlin52.emplace_back(650.0, 1130.0);
        berlin52.emplace_back(1605.0, 620.0);
        berlin52.emplace_back(1220.0, 580.0);
        berlin52.emplace_back(1465.0, 200.0);
        berlin52.emplace_back(1530.0, 5.0);
        berlin52.emplace_back(845.0, 680.0);
        berlin52.emplace_back(725.0, 370.0);
        berlin52.emplace_back(145.0, 665.0);
        berlin52.emplace_back(415.0, 635.0);
        berlin52.emplace_back(510.0, 875.0);
        berlin52.emplace_back(560.0, 365.0);
        berlin52.emplace_back(300.0, 465.0);
        berlin52.emplace_back(520.0, 585.0);
        berlin52.emplace_back(480.0, 415.0);
        berlin52.emplace_back(835.0, 625.0);
        berlin52.emplace_back(975.0, 580.0);
        berlin52.emplace_back(1215.0, 245.0);
        berlin52.emplace_back(1320.0, 315.0);
        berlin52.emplace_back(1250.0, 400.0);
        berlin52.emplace_back(660.0, 180.0);
        berlin52.emplace_back(410.0, 250.0);
        berlin52.emplace_back(420.0, 555.0);
        berlin52.emplace_back(575.0, 665.0);
        berlin52.emplace_back(1150.0, 1160.0);
        berlin52.emplace_back(700.0, 580.0);
        berlin52.emplace_back(685.0, 595.0);
        berlin52.emplace_back(685.0, 610.0);
        berlin52.emplace_back(770.0, 610.0);
        berlin52.emplace_back(795.0, 645.0);
        berlin52.emplace_back(720.0, 635.0);
        berlin52.emplace_back(760.0, 650.0);
        berlin52.emplace_back(475.0, 960.0);
        berlin52.emplace_back(95.0, 260.0);
        berlin52.emplace_back(875.0, 920.0);
        berlin52.emplace_back(700.0, 500.0);
        berlin52.emplace_back(555.0, 815.0);
        berlin52.emplace_back(830.0, 485.0);
        berlin52.emplace_back(1170.0, 65.0);
        berlin52.emplace_back(830.0, 610.0);
        berlin52.emplace_back(605.0, 625.0);
        berlin52.emplace_back(595.0, 360.0);
        berlin52.emplace_back(1340.0, 725.0);
        berlin52.emplace_back(1740.0, 245.0);
    };

    Tour best_tour;
    {
        best_tour.emplace_back(berlin52[1 - 1]);
        best_tour.emplace_back(berlin52[49 - 1]);
        best_tour.emplace_back(berlin52[32 - 1]);
        best_tour.emplace_back(berlin52[45 - 1]);
        best_tour.emplace_back(berlin52[19 - 1]);
        best_tour.emplace_back(berlin52[41 - 1]);
        best_tour.emplace_back(berlin52[8 - 1]);
        best_tour.emplace_back(berlin52[9 - 1]);
        best_tour.emplace_back(berlin52[10 - 1]);
        best_tour.emplace_back(berlin52[43 - 1]);
        best_tour.emplace_back(berlin52[33 - 1]);
        best_tour.emplace_back(berlin52[51 - 1]);
        best_tour.emplace_back(berlin52[11 - 1]);
        best_tour.emplace_back(berlin52[52 - 1]);
        best_tour.emplace_back(berlin52[14 - 1]);
        best_tour.emplace_back(berlin52[13 - 1]);
        best_tour.emplace_back(berlin52[47 - 1]);
        best_tour.emplace_back(berlin52[26 - 1]);
        best_tour.emplace_back(berlin52[27 - 1]);
        best_tour.emplace_back(berlin52[28 - 1]);
        best_tour.emplace_back(berlin52[12 - 1]);
        best_tour.emplace_back(berlin52[25 - 1]);
        best_tour.emplace_back(berlin52[4 - 1]);
        best_tour.emplace_back(berlin52[6 - 1]);
        best_tour.emplace_back(berlin52[15 - 1]);
        best_tour.emplace_back(berlin52[5 - 1]);
        best_tour.emplace_back(berlin52[24 - 1]);
        best_tour.emplace_back(berlin52[48 - 1]);
        best_tour.emplace_back(berlin52[38 - 1]);
        best_tour.emplace_back(berlin52[37 - 1]);
        best_tour.emplace_back(berlin52[40 - 1]);
        best_tour.emplace_back(berlin52[39 - 1]);
        best_tour.emplace_back(berlin52[36 - 1]);
        best_tour.emplace_back(berlin52[35 - 1]);
        best_tour.emplace_back(berlin52[34 - 1]);
        best_tour.emplace_back(berlin52[44 - 1]);
        best_tour.emplace_back(berlin52[46 - 1]);
        best_tour.emplace_back(berlin52[16 - 1]);
        best_tour.emplace_back(berlin52[29 - 1]);
        best_tour.emplace_back(berlin52[50 - 1]);
        best_tour.emplace_back(berlin52[20 - 1]);
        best_tour.emplace_back(berlin52[23 - 1]);
        best_tour.emplace_back(berlin52[30 - 1]);
        best_tour.emplace_back(berlin52[2 - 1]);
        best_tour.emplace_back(berlin52[7 - 1]);
        best_tour.emplace_back(berlin52[42 - 1]);
        best_tour.emplace_back(berlin52[21 - 1]);
        best_tour.emplace_back(berlin52[17 - 1]);
        best_tour.emplace_back(berlin52[3 - 1]);
        best_tour.emplace_back(berlin52[18 - 1]);
        best_tour.emplace_back(berlin52[31 - 1]);
        best_tour.emplace_back(berlin52[22 - 1]);
        best_tour.emplace_back(berlin52[1 - 1]);
    }

    std::cout << "Best tour length: " << k::fitness(best_tour) << std::endl;
    printTour(berlin52, best_tour);
    plot::city_graph(best_tour);

    auto population = createRandomPopulation(berlin52);

    std::cout << "Begin best tour length: " << k::fitness(min(population).tour) << std::endl;
    printTour(berlin52, min(population).tour);
    plot::city_graph(min(population).tour);

    auto result = _main(population, k::P {0.7, 0.01});

    std::cout << "Final best tour length: " << k::fitness(result.tour) << std::endl;
    printTour(berlin52, result.tour);
    plot::city_graph(min(population).tour);
}

}  // namespace lab3

#endif //LAB3_HPP
