#include "Config.hpp"
#include <fstream>
#include <algorithm>
#include <memory>

Config::Config(std::string fileName)
{
    std::ifstream file(fileName);
    std::string line;
    while (std::getline(file, line))
    {
        // FIXME: don't remove # when it's inside a string
        if (line.find('#') != std::string::npos)
            line.erase(line.find('#'), line.size());

        std::string new_str = "";
        bool delete_spaces = true, start = false;
        for (unsigned int i = 0; i < line.size(); ++i)
        {
            if (line[i] == '\"')
            {
                start ? start = false : start = true;
                if (start)
                    delete_spaces = false;
            }
            if (!start)
                delete_spaces = true;
            if (delete_spaces)
            {
                if (line[i] != ' ')
                    new_str += line[i];
            }
            else
                new_str += line[i];
        }
        line = new_str;

        if (line.empty())
            continue;
        std::string key = line.substr(0, line.find('='));
        std::string value = line.substr(line.find('=') + 1, line.size());
        if (value[0] == '"')
        {
            value.erase(value.begin());
            value.erase(value.end() - 1);
            values.emplace(key, value);
        }
        else if (value.find_first_of(".eE") != std::string::npos)
            values.emplace(key, std::stod(value));
        else
            values.emplace(key, (unsigned int) std::stoi(value));
    }
    file.close();
}

unsigned int Config::get_int(std::string key, unsigned int defaultValue)
{
    if (values.find(key) == values.end())
        values.emplace(key, defaultValue);
    return values[key].i;
}

double Config::get_double(std::string key, double defaultValue)
{
    if (values.find(key) == values.end())
        values.emplace(key, defaultValue);
    return values[key].d;
}

std::string Config::get_string(std::string key, std::string defaultValue)
{
    if (values.find(key) == values.end())
        values.emplace(key, defaultValue);
    return *values[key].s;
}