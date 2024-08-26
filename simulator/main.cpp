#include "StarArchitecture.hpp"
#include "Config.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <regex>

std::pair<unsigned int, std::vector<Gate *>> readGatesFromFile(std::string fileName)
{
    std::ifstream file(fileName);
    std::vector<Gate *> gates;
    unsigned int numGates, numQubits = 0;

    file >> numGates;
    for (unsigned int id = 0; id < numGates; ++id)
    {
        std::string gateType;
        file >> gateType;
        if (gateType == "rz")
        {
            unsigned int qubit;
            double angle;
            file >> qubit >> angle;
            if (fabs(fmod(angle, M_PI)) < 1e-6)
                continue;
            numQubits = std::max(numQubits, qubit + 1);
            gates.push_back(new RzGate(id, qubit, angle));
        }
        else if (gateType == "h")
        {
            unsigned int qubit;
            file >> qubit;
            numQubits = std::max(numQubits, qubit + 1);
            gates.push_back(new HadamardGate(id, qubit));
        }
        else if (gateType == "cx")
        {
            unsigned int control, target;
            file >> control >> target;
            numQubits = std::max(numQubits, std::max(control, target) + 1);
            gates.push_back(new CNOTGate(id, control, target));
        }
    }
    return {numQubits, gates};
}

void writeHeatmap(std::string fileName, std::map<unsigned int, std::map<unsigned int, double>> &data, unsigned int numberOfRuns)
{
    std::ofstream file(fileName);
    for (auto [epoch, timeStepHeatmap] : data)
        for (auto [i, heat] : timeStepHeatmap)
            file << heat / numberOfRuns << ",\n"[i == timeStepHeatmap.size() - 1];
    file.close();
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <config file>" << '\n';
        return 1;
    }

    Config config(argv[1]);
    auto inputPath = config.get_string("input_dir");
    auto inputFileRegex = std::regex(config.get_string("input_file", "gates_.*"));
    auto outputsFolder = std::filesystem::path(config.get_string("output_dir"));
    unsigned int codeDistance = config.get_int("code_distance");
    double physicalQubitErrorRate = config.get_double("physical_qubit_error_rate");
    bool debug = config.get_int("debug");
    unsigned int numberOfRuns = config.get_int("number_of_runs", 100);
    auto compiler = config.get_string("compiler", "dynamic");
    double compressionFactor = config.get_double("compression_factor", 0);
    double mstComputationFrequency = config.get_double("mst_computation_frequency", 100);

    auto rotationErrorModel = [](double physicalQubitErrorRate, double) -> double
    {
        return physicalQubitErrorRate;
    };

    std::cout << "Output directory: " << outputsFolder << '\n';
    std::string command = "mkdir -p \"" + outputsFolder.string() + "\"";
    system(command.c_str());
    command = "cp " + std::string(argv[1]) + " \"" + outputsFolder.string() + "/config.cfg\"";
    system(command.c_str());

    for (auto &file : std::filesystem::directory_iterator(inputPath))
    {
        if (std::regex_match(file.path().filename().string(), inputFileRegex) == false)
            continue;
        auto outputFolder = outputsFolder / file.path().filename();
        std::string command = "mkdir -p \"" + outputFolder.string() + "\"";
        system(command.c_str());
        std::ofstream logFile(outputFolder / std::filesystem::path("log"));
        logFile << file.path() << '\n';
        std::cout << file.path() << '\n';

        auto [numQubits, gates] = readGatesFromFile(file.path());
        unsigned int numQubitsPerRow = std::ceil(std::sqrt(numQubits));
        unsigned int numQubitsPerColumn = std::ceil((double)numQubits / numQubitsPerRow);
        numQubits = numQubitsPerRow * numQubitsPerColumn;

        std::map<unsigned int, std::map<unsigned int, double>> dataQubitHeatmap;
        std::map<unsigned int, std::map<unsigned int, double>> ancillaQubitHeatmap;
        BaseStarArchitecture::initStaticStuff(debug, logFile);
        double averageTimeTaken = 0;
        for (unsigned int seed = 0; seed < numberOfRuns; seed++)
        {
            BaseStarArchitecture *starArchitecture = nullptr;
            if (compiler == "static")
                starArchitecture = new StaticStarArchitecture(numQubitsPerRow * 2, numQubitsPerColumn * 2, numQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, compressionFactor);
            else if (compiler == "dynamic")
                starArchitecture = new DynamicStarArchitecture(numQubitsPerRow * 2, numQubitsPerColumn * 2, numQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, mstComputationFrequency, compressionFactor);
            else if (compiler == "autobraid")
                starArchitecture = new AutoBraidStarArchitecture(numQubitsPerRow * 2, numQubitsPerColumn * 2, numQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, compressionFactor);
            else
            {
                std::cout << "Invalid compiler: " << compiler << '\n';
                return 1;
            }
            starArchitecture->addGates(gates);

            starArchitecture->simulate();
            starArchitecture->updateQubitsBusyTimes(dataQubitHeatmap, false);
            starArchitecture->updateQubitsBusyTimes(ancillaQubitHeatmap, true);
            logFile << "Done in " << starArchitecture->getGlobalTime() << " cycles" << std::endl;
            averageTimeTaken += starArchitecture->getGlobalTime();
        }
        averageTimeTaken /= numberOfRuns;

        writeHeatmap(outputFolder / std::filesystem::path("dataq_heatmap.csv"), dataQubitHeatmap, numberOfRuns);
        writeHeatmap(outputFolder / std::filesystem::path("ancillaq_heatmap.csv"), ancillaQubitHeatmap, numberOfRuns);

        outputFolder = outputFolder / std::filesystem::path("logs");
        command = "mkdir -p \"" + outputFolder.string() + "\"";
        system(command.c_str());
        BaseStarArchitecture::logTimeTaken(outputFolder);
        logFile << "Total execution average: " << averageTimeTaken << '\n';
        BaseStarArchitecture::resetStaticStuff();
    }

    return 0;
}