#ifndef RANDOM_HPP
#define RANDOM_HPP
#include <random>

namespace random
{
int _int(int l, int r)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(l, r);
    return dist(gen);
}
float _float(float l, float r)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(l, r);
    return (float) dist(gen);
}
float _double(double l, double r)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(l, r);
    return dist(gen);
}
}

#endif //RANDOM_HPP
