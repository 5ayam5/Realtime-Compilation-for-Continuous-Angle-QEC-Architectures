#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <unordered_map>
#include <memory>

class Config
{
public:
    Config(std::string fileName);

    unsigned int get_int(std::string key, unsigned int defaultValue = 0);
    double get_double(std::string key, double defaultValue = 0.0);
    std::string get_string(std::string key, std::string defaultValue = "");

private:
    struct Value
    {
        unsigned int i;
        double d;
        std::unique_ptr<std::string> s;

        Value() {};
        Value(unsigned int i) : i(i), d(i) {};
        Value(double d) : d(d) {};
        Value(std::string s) : s(std::make_unique<std::string>(s)) {};
        ~Value() {};
    };

    std::unordered_map<std::string, Value> values;
};

#endif