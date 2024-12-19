#ifndef LAB4_HPP
#define LAB4_HPP
#include <utility>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <memory>
#include <iostream>
#include "Random.hpp"
#include "Bocchiemark.hpp"

namespace lab4
{

namespace k
{
std::unordered_set<std::string> term_set
    {
        "x1",
        "x2",
        "x3",
        "x4",
        "x5",
        "4000",
        "1"
    };
std::unordered_set<std::string> func_set
    {
        "+",
        "-",
        "*",
        "/",
        "abs",
        "sin",
        "cos",
        "exp",
        "pow"
    };

const uint16_t MAX_DEPTH = 200;
const uint16_t MAX_UNARY_DEPTH = 10;
const uint16_t MAX_ATTEMPTS = 10;
}  // namespace k

struct Node
{
    enum class Type { kTerm, kFunc } type;
    enum class FuncType { kBinary, kUnary } func_type;
    std::string value;

    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;

    Node(Type type, std::string value)
        : type(type)
        , value(std::move(value))
        , left(nullptr)
        , right(nullptr)
    {
        if (this->value == "+") func_type = FuncType::kBinary;
        if (this->value == "-") func_type = FuncType::kBinary;
        if (this->value == "*") func_type = FuncType::kBinary;
        if (this->value == "/") func_type = FuncType::kBinary;
        if (this->value == "abs") func_type = FuncType::kUnary;
        if (this->value == "sin") func_type = FuncType::kUnary;
        if (this->value == "cos") func_type = FuncType::kUnary;
        if (this->value == "exp") func_type = FuncType::kUnary;
        if (this->value == "pow") func_type = FuncType::kBinary;
    }

    bool isBinary()
    {
        return this->func_type == FuncType::kBinary;
    }

    // Собираем все узлы в вектор
    static void collectNodes(const std::shared_ptr<Node>& node, std::vector<std::shared_ptr<Node>>& nodes) {
        if (!node) return;

        nodes.push_back(node);

        if (node->type == Node::Type::kFunc) {
            collectNodes(node->left, nodes);
            if (node->func_type == Node::FuncType::kBinary)
                collectNodes(node->right, nodes);
        }
    }

    // Функция для получения случайного узла в дереве
    static std::shared_ptr<Node> getRandomNode(const std::shared_ptr<Node>& root) {
        std::vector<std::shared_ptr<Node>> nodes;
        collectNodes(root, nodes);
        return nodes[random::_int(0, nodes.size() - 1)];
    }

    // Функция для вычисления глубины дерева
    static int calculateDepth(const std::shared_ptr<Node>& root) {
        if (!root) return 0;
        if (root->type == Node::Type::kTerm) return 1;
        int left_depth = calculateDepth(root->left);
        int right_depth = (root->func_type == Node::FuncType::kBinary) ? calculateDepth(root->right) : 0;
        return std::max(left_depth, right_depth) + 1;
    }

    static std::shared_ptr<Node> generateTreeWithTerms(
        std::unordered_set<std::shared_ptr<Node>>& terms,
        int maxUnaryDepth,
        int currentUnaryDepth = 0)
    {
        if (terms.size() == 1) {
            // Если остался только один термин, он становится листом
            auto term = *terms.begin();
            terms.erase(terms.begin());
            return term;
        }

        // Выбираем случайный оператор
        auto func = *std::next(k::func_set.begin(), random::_int(0, k::func_set.size() - 1));
        std::shared_ptr<Node> node = std::make_shared<Node>(Node::Type::kFunc, func);

        if (node->isBinary()) {
            // Если оператор бинарный, разделяем термы на две группы
            size_t split = terms.size() / 2;
            std::unordered_set<std::shared_ptr<Node>> leftTerms, rightTerms;

            auto it = terms.begin();
            for (size_t i = 0; i < split; ++i, ++it) {
                leftTerms.insert(*it);
            }
            for (; it != terms.end(); ++it) {
                rightTerms.insert(*it);
            }

            // Рекурсивно создаём поддеревья
            node->left = generateTreeWithTerms(leftTerms, maxUnaryDepth);
            node->right = generateTreeWithTerms(rightTerms, maxUnaryDepth);
        } else {
            // Если оператор унарный, проверяем ограничение по глубине унарных операторов
            if (currentUnaryDepth >= maxUnaryDepth) {
                // Превышение глубины унарных операторов: выбираем термин
                auto term = *terms.begin();
                terms.erase(terms.begin());
                return term;
            }

            // Рекурсивно создаём поддерево для унарного оператора
            node->left = generateTreeWithTerms(terms, maxUnaryDepth, currentUnaryDepth + 1);
        }

        return node;
    }

    static std::shared_ptr<Node> generateCompleteTreeWithTerms(std::unordered_set<std::string> term_set = k::term_set, int depth = k::MAX_DEPTH, int maxUnaryDepth = k::MAX_UNARY_DEPTH) {
        // Собираем все элементы из term_set в terms
        std::unordered_set<std::shared_ptr<Node>> terms;
        for (const auto& term : term_set) {
            terms.insert(std::make_shared<Node>(Node::Type::kTerm, term));
        }

        // Генерируем дерево с гарантией, что все термины будут включены
        return generateTreeWithTerms(terms, maxUnaryDepth);
    }

    static std::unordered_set<std::string> collectTerms(const std::shared_ptr<Node>& root) {
        std::unordered_set<std::string> collected_terms;

        // Рекурсивный обход поддерева
        std::function<void(const std::shared_ptr<Node>&)> traverse = [&](const std::shared_ptr<Node>& node) {
            if (!node) return;

            if (node->type == Node::Type::kTerm && k::term_set.count(node->value)) {
                collected_terms.insert(node->value);
            }

            if (node->left) traverse(node->left);
            if (node->right) traverse(node->right);
        };

        traverse(root);
        return collected_terms;
    }

};

struct Entity
{
private:
    // Рекурсивное вычисление значения дерева
    static double evaluate(const std::shared_ptr<Node>& node, const std::vector<double>& x) {
        if (! node) return 0xFFFFFFFF;
        try {
            if (node->type == Node::Type::kTerm) {
                // Если это терминальная переменная (например x1), извлекаем значение из вектора
                if (node->value[0] == 'x') {
                    // Извлекаем индекс переменной (например "x1" -> 1)
                    int index = std::stoi(node->value.substr(1)) - 1; // x1 -> индекс 0
                    if (index >= 0 && index < x.size()) {
                        return x[index];
                    }
                    std::cerr << "[Error] Variable index out of range: " << node->value << "\n";
                    return 0xFFFFFFFF;
                } else {
                    // Если это число, просто возвращаем его значение
                    try {
                        return std::stod(node->value); // Преобразуем строку в число
                    } catch (const std::invalid_argument& e) {
                        std::cerr << "[Error] Invalid number format in node value: " << node->value << "\n";
                        return 0xFFFFFFFF;
                    }
                }
            } else if (node->type == Node::Type::kFunc) {
                if (node->func_type == Node::FuncType::kUnary)
                    return applyUnaryOperator(node->value, evaluate(node->left, x));
                else if (node->func_type == Node::FuncType::kBinary)
                    return applyBinaryOperator(node->value, evaluate(node->left, x), evaluate(node->right, x));
            }
            std::cerr << "[Error] Invalid node type encountered!\n";
        } catch (const std::exception& e) {
            std::cerr << "[Error] Exception caught: " << e.what() << "\n";
        } catch (...) {
            std::cerr << "[Error] Unknown exception caught!\n";
        }

        // Возврат значения ошибки в случае любой ошибки
        return 0xFFFFFFFF;
    }

    // Применение бинарного оператора
    static double applyBinaryOperator(const std::string& op, double left, double right) {
        if (op == "+") return left + right;
        if (op == "-") return left - right;
        if (op == "*") return left * right;
        if (op == "/") return right != 0 ? left / right : 1e9; // Избегаем деления на 0
        if (op == "pow") return std::pow(left, right);
        throw std::invalid_argument("Unsupported binary operator: " + op);
    }

    // Применение унарного оператора
    static double applyUnaryOperator(const std::string& op, double operand) {
        if (op == "abs") return std::abs(operand);
        if (op == "sin") return std::sin(operand);
        if (op == "cos") return std::cos(operand);
        if (op == "exp") return std::exp(operand);
        throw std::invalid_argument("Unsupported unary operator: " + op);
    }

public:
    std::shared_ptr<Node> root;

    Entity() : root(nullptr) {}
    Entity(std::shared_ptr<Node> node) : root(std::move(node)) {}

    double operator()(const std::vector<double>& x) const {
        if (!root) throw std::runtime_error("Tree is empty!");
        return evaluate(root, x);
    }
};


typedef std::vector<Entity> Population;
typedef std::vector<double> Fitness;

namespace k
{
const int16_t MIN = -600;
const int16_t MAX = 600;

const uint16_t F = 1'000;

struct P
{
    float _CROSSING = 0.5;     // Вероятность скрещивания
    float _MUTATION = 0.0001;  // Вероятность мутации
};

const uint16_t N = 100; // размер популяции
const uint8_t n = 5;

// f(x) = sum( x(i) ^ 2 / 4000) - prod(cos(x(i)) / sqrt(i))) + 1
// i = 1:5
// -600 <= x(i) <= 600
double f(const std::vector<double>& x) {

    if (x.size() != 5) {
        throw std::invalid_argument("Input vector must have exactly 5 elements.");
    }

    double sum = 0.0;
    double prod = 1.0;

    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] < -600 || x[i] > 600) {
            throw std::out_of_range("Each element of x must be in the range [-600, 600].");
        }

        double xi = x[i];
        size_t index = i + 1; // Индексы начинаются с 1 по условию задачи

        sum += (xi * xi) / 4000.0;
        prod *= std::cos(xi / std::sqrt(static_cast<double>(index)));
    }

    return sum - prod + 1.0;
}

// модуль разности
double fitness(const Entity& entity)
{
    static auto x_list = [](){
        std::vector<std::vector<double>> x_list;
        for (size_t i = 0; i < F; i++) {
            std::vector<double> x;
            for (size_t j = 0; j < n; ++j)
                x.push_back(random::_int(MIN, MAX));
            x_list.push_back(x);
        }
        return x_list;
    }();

    double sum = 0.0;
    for (const auto& x : x_list)
        sum += std::abs(entity(x) - f(x));
    return sum;
}

Fitness fitness(const Population& population)
{
    Fitness out;

    for (const auto& entity: population)
        out.push_back(fitness(entity));

    return out;
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
// Кроссинговер поддеревьев
Population crossover(const Population &population) {
    Population offspring;

    for (size_t i = 0; i < population.size() / 2; i++) {
        const auto& parent1 = population[i];
        const auto& parent2 = population[i + population.size() / 2];

        bool crossover_performed = false;

        // Попытки найти подходящие узлы
        for (int attempt = 0; attempt < k::MAX_ATTEMPTS; ++attempt) {
            std::shared_ptr<Node> node1 = Node::getRandomNode(parent1.root);
            std::shared_ptr<Node> node2 = Node::getRandomNode(parent2.root);

            if (node1->type == node2->type) {
                // Обмениваемся поддеревьями
                if (node1->type == Node::Type::kFunc) {
                    std::swap(node1->left, node2->left);
                    std::swap(node1->right, node2->right);
                } else if (node1->type == Node::Type::kTerm) {
                    std::swap(node1->value, node2->value);
                }

                // Проверяем глубину потомков
                int depth1 = Node::calculateDepth(parent1.root);
                int depth2 = Node::calculateDepth(parent2.root);

                if (depth1 <= k::MAX_DEPTH) {
                    offspring.push_back(parent1);
                } else {
                    offspring.push_back(parent2);
                }

                if (depth2 <= k::MAX_DEPTH) {
                    offspring.push_back(parent2);
                } else {
                    offspring.push_back(parent1);
                }

                crossover_performed = true;
                break;
            }
        }

        // Если подходящие узлы так и не найдены, добавляем родителей без изменений
        if (!crossover_performed) {
            offspring.push_back(parent1);
            offspring.push_back(parent2);
        }
    }

    return offspring;
}



// Растущая мутация
Population mutate(const Population &population, float mutation_rate) {
    Population mutated_population = population;

    for (auto& entity : mutated_population) {
        if (random::_int(0, 100) < mutation_rate * 100) {  // С вероятностью mutation_rate
            std::shared_ptr<Node> node_to_mutate = Node::getRandomNode(entity.root);

            if (node_to_mutate->type == Node::Type::kFunc) {
                // Сбор терминов из поддерева, которое будет заменено
                std::unordered_set<std::string> subtree_terms = Node::collectTerms(node_to_mutate);

                // Генерируем новое случайное поддерево
                node_to_mutate->left = nullptr;
                if (node_to_mutate->func_type == Node::FuncType::kBinary)
                    node_to_mutate->right = nullptr;

                node_to_mutate->left = Node::generateCompleteTreeWithTerms(subtree_terms);
                if (node_to_mutate->func_type == Node::FuncType::kBinary)
                    node_to_mutate->right = Node::generateCompleteTreeWithTerms(subtree_terms);
            } else {
                // Если узел терминальный, заменяем его другим терминальным узлом
                node_to_mutate->value = *std::next(k::term_set.begin(), random::_int(0, k::term_set.size() - 1));
            }

            // Проверяем глубину дерева после мутации
            if (Node::calculateDepth(entity.root) > k::MAX_DEPTH) {
                entity.root = Node::generateCompleteTreeWithTerms();  // Если превышает, создаем новое дерево
            }
        }
    }

    return mutated_population;
}

}  // namespace op

// Функция для печати дерева
void drawTree(const std::shared_ptr<Node>& node, const std::string& prefix = "", bool isLeft = true) {
    if (!node) return;

    std::cout << prefix;

    if (isLeft) {
        std::cout << "|-- ";
    } else {
        std::cout << "`-- ";
    }

    // Выводим сам узел
    std::cout << node->value << std::endl;

    // Рекурсивно выводим поддеревья
    if (node->type == Node::Type::kFunc) {
        if (node->func_type == Node::FuncType::kUnary) {
            drawTree(node->left, prefix + (isLeft ? "|   " : "    "), true);
        } else if (node->func_type == Node::FuncType::kBinary) {
            drawTree(node->left, prefix + (isLeft ? "|   " : "    "), true);
            drawTree(node->right, prefix + (isLeft ? "|   " : "    "), false);
        }
    }
}

auto min = [](const Population& population)
{
    return *std::min_element(population.begin(), population.end(),
                             [](const auto &a, const auto &b)
                             {
                                 return k::fitness(a) < k::fitness(b);
                             });
};

Entity _main(Population population, k::P p)
{
    size_t epoch_counter = 0;

    for (size_t i = 0; i < 1'000; ++i)
    {
        epoch_counter++;
        try {

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

            if (i % 10 == 0)
            {
                std::cout << "Epoch " << i << ", avg fitness (abs diff): " << k::fitness(min(population)) / k::F << std::endl;
            }
            if (i % 100 == 0)
            {
                drawTree(min(population).root);
            }
            }
        catch (const std::exception& e)
        {
            std::cout << "Aboba: " << e.what() << std::endl;
        }
        catch (...)
        {
            std:: cout << "Somthing went wrong ...\n" ;
        };
    }

    std::cout << "Educated by " << epoch_counter << " epochs\n";

    return min(population);
}

void main()
{

    while (false) {
        Entity _;
        _.root = Node::generateCompleteTreeWithTerms();

        drawTree(_.root);

        std::cout << "Aboba: " << _({1,2,3,4,5}) << std::endl;
        std::cout << "Fitness: " << k::fitness(_) << std::endl;

        return;
    }

    bocchie::mark mark("run");

    Population population;
    for (size_t i = 0; i < k::N; ++i)
        population.emplace_back(Node::generateCompleteTreeWithTerms());

    auto best = mark.run(_main, population, k::P {0.6, 0.1});

    drawTree(best.root);

    std::cout << mark.to_json<bocchie::accuracy::milliseconds>() << std::endl;

}


}  // namespace lab4

#endif //LAB4_HPP
