#ifndef BOCCHIEMARK_BOCCHIEMARK_H
#define BOCCHIEMARK_BOCCHIEMARK_H
#include <string>
#include <cstdint>
#include <sstream>
#include <chrono>
#include <functional>
#include <map>

namespace bocchie
{
class mark;
namespace instances
{
std::map<std::string /* <- runnable_name */, mark> map;
}

enum class accuracy {
    seconds,
    milliseconds,
    microseconds,
    nanoseconds
};

template<bocchie::accuracy accuracy>
struct SillyTime;

template<>
struct SillyTime<bocchie::accuracy::seconds> {
    static double cast(const uint64_t& time_nanosec) noexcept {
        return static_cast<double>(time_nanosec) / 1e9;
    }
    static std::string to_string() noexcept {
        return "seconds";
    }
};

template<>
struct SillyTime<bocchie::accuracy::milliseconds> {
    static double cast(const uint64_t& time_nanosec) noexcept {
        return static_cast<double>(time_nanosec) / 1e6;
    }
    static std::string to_string() noexcept {
        return "milliseconds";
    }
};

template<>
struct SillyTime<bocchie::accuracy::microseconds> {
    static double cast(const uint64_t& time_nanosec) noexcept {
        return static_cast<double>(time_nanosec) / 1e3;
    }
    static std::string to_string() noexcept {
        return "microseconds";
    }
};

template<>
struct SillyTime<bocchie::accuracy::nanoseconds> {
    static uint64_t cast(const uint64_t& time_nanosec) noexcept {
        return time_nanosec;
    }
    static std::string to_string() noexcept {
        return "nanoseconds";
    }
};

class mark {
public:
    friend class mark_test; // lol, I think I`ll implement it hahaha
private:
    std::string runnable_name;
    uint64_t last_time_nanosec;
    uint64_t min_time_nanosec;
    uint64_t max_time_nanosec;
    uint64_t total_time_nanosec;
    uint64_t total_runs;
public:
    explicit mark(std::string runnable_name) :
        runnable_name(std::move(runnable_name)),
        last_time_nanosec(0x00000000),
        min_time_nanosec(UINT64_MAX),
        max_time_nanosec(0x00000000),
        total_time_nanosec(0x00000000),
        total_runs(0x00000000)
    {}

    explicit mark():
        last_time_nanosec(0x00000000),
        min_time_nanosec(UINT64_MAX),
        max_time_nanosec(0x00000000),
        total_time_nanosec(0x00000000),
        total_runs(0x00000000)
    {}

    template<class... Args>
    void run(void (*func)(Args...), Args... args)
    {
        auto start = std::chrono::high_resolution_clock::now();
        try {
            func(args...);
        } catch (std::exception& e) {
            throw (std::string("Running \"" + this->get_runnable_name() + "\" threw exception: " + e.what() + "\n"));
        }
        auto end = std::chrono::high_resolution_clock::now();
        this->last_time_nanosec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();\
if (this->last_time_nanosec < this->min_time_nanosec) this->min_time_nanosec = this->last_time_nanosec;
        if (this->last_time_nanosec > this->max_time_nanosec) this->max_time_nanosec = this->last_time_nanosec;
        this->total_time_nanosec += this->last_time_nanosec;
        this->total_runs++;
    }

    template<class Func, class... Args>
    auto run(Func func, Args... args)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = func(args...); // FIXME: mb call in try-catch
        auto end = std::chrono::high_resolution_clock::now();
        this->last_time_nanosec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();\
if (this->last_time_nanosec < this->min_time_nanosec) this->min_time_nanosec = this->last_time_nanosec;
        if (this->last_time_nanosec > this->max_time_nanosec) this->max_time_nanosec = this->last_time_nanosec;
        this->total_time_nanosec += this->last_time_nanosec;
        this->total_runs++;
        return static_cast<decltype(result)>(result);
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] std::string to_json() const
    {
        auto accuracy_string = SillyTime<accuracy>::to_string();
        std::stringstream ss;
        ss << "{\n"
           << "\t\"runnable_name\": \"" << this->runnable_name << "\",\n"
           << "\t\"min_time_" << accuracy_string << "\": " << this->get_min_time<accuracy>() << ",\n"
           << "\t\"max_time_" << accuracy_string << "\": " << this->get_max_time<accuracy>() << ",\n"
           << "\t\"last_time_" << accuracy_string << "\": " << this->get_last_time<accuracy>() << ",\n"
           << "\t\"total_time_" << accuracy_string << "\": " << this->get_total_time<accuracy>() << ",\n"
           << "\t\"total_runs" << "\": " << this->get_total_runs() << ",\n"
           << "\t\"avg_time_" << accuracy_string << "\": " << this->get_avg_time<accuracy>() << "\n"
           << "}";
        return ss.str();
    }

    [[nodiscard]] std::string get_runnable_name() const noexcept
    {
        return this->runnable_name;
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] auto get_last_time() const noexcept
    {
        return SillyTime<accuracy>::cast(this->last_time_nanosec);
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] auto get_min_time() const noexcept
    {
        return SillyTime<accuracy>::cast(this->min_time_nanosec);
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] auto get_max_time() const noexcept
    {
        return SillyTime<accuracy>::cast(max_time_nanosec);
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] auto get_total_time() const noexcept
    {
        return SillyTime<accuracy>::cast(total_time_nanosec);
    }

    [[nodiscard]] uint64_t get_total_runs() const noexcept
    {
        return this->total_runs;
    }

    template<bocchie::accuracy accuracy>
    [[nodiscard]] auto get_avg_time() const noexcept
    {
        return total_runs == 0 ? this->get_total_time<accuracy>() :
               SillyTime<accuracy>::cast(total_time_nanosec / total_runs);
    }
};

class manager
{
public:
    static mark& get(const std::string& name)
    {
        return instances::map[name];
    }
};
}
#endif //BOCCHIEMARK_BOCCHIEMARK_H