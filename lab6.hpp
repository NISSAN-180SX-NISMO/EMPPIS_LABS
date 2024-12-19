#ifndef LAB6_HPP
#define LAB6_HPP

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <set>

namespace lab6
{

typedef std::pair<float, float> City;
typedef std::vector<City> Tour;
typedef uint16_t CityNumber;
typedef std::map<CityNumber, City> CityMap;
typedef std::map<City, CityNumber> InvCityMap;

using CityPair = std::pair<CityNumber, CityNumber>; // Пара индексов городов
typedef std::map<CityPair, double> CityGraph; // Мапа для хранения расстояний / флогистона между городами


struct Ant {
    Tour tour;
    float tourLength;
};

typedef std::vector<Ant> Population;

namespace k
{

const size_t N = 1000;
const float ALPHA = 1.0;  // Влияние феромона
const float BETA = 4.0;   // Влияние эвристики
const float RHO = 0.6;    // Испарение феромона
const float Q = 1.0;    // Количество феромона, оставляемого муравьем

float calculateDistance(const City& a, const City& b) {
    float dx = a.first - b.first;
    float dy = a.second - b.second;
    return std::sqrt(dx * dx + dy * dy);
}

float fitness(const Tour& tour) {
    if (tour.size() < 2) return 0.0f;

    float totalDistance = 0.0f;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        const City& from = tour[i];
        const City& to = tour[i + 1];
        totalDistance += calculateDistance(from, to);
    }
    return totalDistance;
}

float fitness(const Ant& ant)
{
    return fitness(ant.tour);
}

std::vector<float> fitness(const Population& population)
{
    std::vector<float> out;
    for (const auto& ant: population)
        out.push_back(fitness(ant));
    return out;
}

}  // namespace k

class AntColonyOptimization {
public:
    AntColonyOptimization(CityMap cityMap,
                          InvCityMap cityInvMap,
                          CityGraph pheromoneGraph,
                          CityGraph heuristicGraph,
                          CityGraph distanceGraph,
                          size_t numAnts, size_t maxIterations)
        : cityMap(std::move(cityMap)),
          cityInvMap(std::move(cityInvMap)),
          pheromoneGraph(std::move(pheromoneGraph)),
          heuristicGraph(std::move(heuristicGraph)),
          distanceGraph(std::move(distanceGraph)),
          numAnts(numAnts),
          maxIterations(maxIterations),
          bestTourLength(std::numeric_limits<float>::max()) {
        ants.resize(numAnts);
    }

    void run() {
        const size_t numCities = cityMap.size();
        std::random_device rd;
        std::mt19937 gen(rd());

        for (size_t iteration = 0; iteration < maxIterations; ++iteration) {

            placeAnts(gen);
            constructTours(gen);
            updatePheromones();
            if (iteration == 0)
            {
                auto best = *std::min_element(this->ants.begin(), this->ants.end(), [](const Ant& a, const Ant& b)
                {
                    return k::fitness(a) < k::fitness(b);
                });
                plot::city_graph(best.tour);
            }
            if (iteration % 10 == 0)
                printIterationResult(iteration);
        }

        printFinalResult();
    }

private:
    CityMap cityMap;
    InvCityMap cityInvMap;
    CityGraph pheromoneGraph;
    CityGraph heuristicGraph;
    CityGraph distanceGraph;
    size_t numAnts;
    size_t maxIterations;
    Population ants;
    Tour bestTour;
    float bestTourLength;

    void placeAnts(std::mt19937& gen) {
        const size_t numCities = cityMap.size();
        for (auto& ant : ants) {
            ant.tour.clear();
            ant.tourLength = 0.0f;
            ant.tour.push_back(cityMap.at(1));
        }
    }

    void constructTours(std::mt19937& gen) {
        for (auto& ant : ants) {
            std::set<CityNumber> visited;
            visited.insert(ant.tour.front().first);

            while (ant.tour.size() < cityMap.size()) {
                CityNumber currentCity = cityInvMap.at(ant.tour.back());
                chooseNextCity(ant, visited, currentCity, gen);
                localPheromoneUpdate(ant); // Локальное обновление феромонов
            }
            ant.tour.push_back(cityMap.at(1));
            ant.tourLength = k::fitness(ant);
            updateBestTour(ant);
        }
    }

    void chooseNextCity(Ant& ant, std::set<CityNumber>& visited, CityNumber currentCity, std::mt19937& gen) {
        std::vector<std::pair<CityNumber, float>> probabilities;
        double probabilitySum = 0.0f;

        for (const auto& [nextCity, _] : cityMap) {
            if (visited.find(nextCity) != visited.end()) continue;

            float pheromone = pheromoneGraph.at({currentCity, nextCity});
            float heuristic = heuristicGraph.at({currentCity, nextCity});
            double probability = std::pow(pheromone, k::ALPHA) * std::pow(heuristic, k::BETA);

            probabilities.emplace_back(nextCity, probability);
            probabilitySum += probability;
        }

        float randomValue = random::_double(0, probabilitySum);

        CityNumber nextCity = 0;
        for (const auto& [city, probability] : probabilities) {
            randomValue -= probability;
            if (randomValue <= 0) {
                nextCity = city;
                break;
            }
        }

        visited.insert(nextCity);
        ant.tour.push_back(cityMap.at(nextCity));
    }

    void localPheromoneUpdate(const Ant& ant) {
        const float localRho = 0.1; // Локальный коэффициент испарения
        for (size_t i = 0; i < ant.tour.size() - 1; ++i) {
            CityNumber from = ant.tour[i].first;
            CityNumber to = ant.tour[i + 1].first;
            pheromoneGraph[{from, to}] *= (1 - localRho);
            pheromoneGraph[{to, from}] *= (1 - localRho);
        }
    }

    void updateBestTour(const Ant& ant) {
        if (ant.tourLength < bestTourLength) {
            bestTour = ant.tour;
            bestTourLength = ant.tourLength;
        }
    }

    void updatePheromones() {
        for (auto& [edge, pheromone] : pheromoneGraph) {
            pheromone *= (k::RHO);
        }

        for (const auto& ant : ants) {
            float pheromoneDeposit = k::Q / ant.tourLength;

            for (size_t i = 0; i < ant.tour.size() - 1; ++i) {
                CityNumber from = ant.tour[i].first;
                CityNumber to = ant.tour[i + 1].first;
                pheromoneGraph[{from, to}] += pheromoneDeposit;
                pheromoneGraph[{to, from}] += pheromoneDeposit;
            }
        }
    }

    void printIterationResult(size_t iteration) const {
        std::cout << "Iteration " << iteration << ": Best length so far = " << bestTourLength << std::endl;
    }

    void printFinalResult() const {
        std::cout << "Final best tour length: " << bestTourLength << std::endl;
        std::cout << "END" << std::endl;
        plot::city_graph(bestTour);
    }
};


int main() {
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

    CityMap berlin52_map = [] (const std::vector<City>& _)
    {
        CityMap __;
        for (size_t number = 1; number <= _.size(); ++number)
            __[number] = _[number - 1];
        return __;
    }(berlin52);

    InvCityMap berlin52_invmap = [] (const std::vector<City>& _)
    {
        InvCityMap __;
        for (size_t number = 1; number <= _.size(); ++number)
            __[_[number - 1]] = number;
        return __;
    }(berlin52);

    CityGraph distanceGraph = [](const CityMap& _)
    {
        CityGraph __;
        for (const auto& col: _)
            for (const auto& line: _)
                __[{col.first, line.first}] = k::calculateDistance(col.second, line.second);
        return __;
    }(berlin52_map);

    CityGraph pheromoneGraph = [&distanceGraph](const CityMap& _)
    {
        CityGraph __;
        for (const auto& col: _)
            for (const auto& line: _)
                __[{col.first, line.first}] = random::_float(0.001, 0.01);
        return __;
    }(berlin52_map);

    CityGraph heuristicGraph = [&distanceGraph](const CityMap& _)
    {
        CityGraph __;
        for (const auto& col: _)
            for (const auto& line: _)
            {
                auto dist = distanceGraph[{col.first, line.first}];
                __[{col.first, line.first}] = dist == 0 ? 0 : 1 / distanceGraph[{col.first, line.first}];
            }
        return __;
    }(berlin52_map);


    AntColonyOptimization aco(berlin52_map, berlin52_invmap, pheromoneGraph, heuristicGraph, distanceGraph, k::N, 100);
    aco.run();


    return 0;
}


}

#endif //LAB6_HPP
