#ifndef __BASESTARARCHITECTURE_HPP__
#define __BASESTARARCHITECTURE_HPP__

#include "Gate.hpp"
#include "LogicalQubit.hpp"
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <functional>
#include <random>
#include <filesystem>
#include <iostream>
#include <boost/functional/hash.hpp>

std::vector<std::pair<unsigned int, unsigned int>> getNeighbourCoordinates(unsigned int row, unsigned int col);

/**
 * @brief Star architecture class as proposed in https://arxiv.org/abs/2303.13181
 *
 */
class BaseStarArchitecture
{
public:
    BaseStarArchitecture() = delete;

    /**
     * @brief get the number of rows in the architecture
     *
     * @return unsigned int
     */
    unsigned int getNumRows() const;

    /**
     * @brief get the number of columns in the architecture
     *
     * @return unsigned int
     */
    unsigned int getNumColumns() const;

    /**
     * @brief get the surface code distance
     *
     * @return unsigned int
     */
    unsigned int getCodeDistance() const;

    /**
     * @brief get the physical qubit error rate
     *
     * @return double
     */
    double getPhysicalQubitErrorRate() const;

    /**
     * @brief get the global time
     *
     * @return double
     */
    double getGlobalTime() const;

    /**
     * @brief get the gates
     *
     * @return const std::vector<Gate>&
     */
    const std::vector<Gate *> &getGates() const;

    /**
     * @brief add a list of gates
     *
     * @param gates the list of gates to add
     */
    void addGates(const std::vector<Gate *> &gates);

    /**
     * @brief add a single gate
     *
     * @param gate the gate to add
     */
    void addGate(Gate *gate);

    /**
     * @brief simulate the architecture until the next gate
     *
     * @return true if there are no more gates to execute
     */
    virtual bool step() = 0;

    /**
     * @brief simulate the architecture until the end
     *
     */
    void simulate();

    /**
     * @brief log (cerr) the _timeTaken map
     *
     */
    static void logTimeTaken(std::filesystem::path outputFolder);

    /**
     * @brief Get the busy times for the qubits between the interval [startTime, endTime]
     *
     * @param startTime the start time of the interval
     * @param endTime the end time of the interval
     * @param ancillas is it for ancillas or data qubits
     * @return std::vector<double> the list of busy times
     */
    void updateQubitsBusyTimes(std::map<unsigned int, std::map<unsigned int, double>> &qubitsBusyTimes, bool ancillas);

    static void initStaticStuff(bool debug, std::ostream &debugStream);

    /**
     * @brief reset the _timeTaken map
     *
     */
    static void resetStaticStuff();

    /**
     * @brief Destroy the BaseStarArchitecture object
     *
     */
    ~BaseStarArchitecture();

protected:
    /**
     * @brief Construct a new BaseStarArchitecture object
     *
     * @param numRows the number of rows in the architecture
     * @param numColumns the number of columns in the architecture
     * @param codeDistance the surface code distance
     * @param physicalQubitErrorRate the physical qubit error rate
     * @param gates the list of gates
     */
    BaseStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor);

    /**
     * @brief initialize a data qubit
     *
     * @param qubitId the id of the data qubit
     * @return the initialized data qubit
     */
    DataQubit *_getDataQubit(unsigned int qubitId);

    /**
     * @brief get the AncillaQubit at row, column
     *
     * @param row the row of the AncillaQubit
     * @param column the column of the AncillaQubit
     * @return AncillaQubit* or nullptr if there is a DataQubit at row, column (or if row, column is out of bounds)
     */
    AncillaQubit *_getAncillaQubit(unsigned int row, unsigned int column);

    /**
     * @brief get the LogicalQubit at row, column
     *
     * @param row the row of the LogicalQubit
     * @param column the column of the LogicalQubit
     * @return LogicalQubit* the LogicalQubit at row, column (requires that either _getDataQubit or _getAncillaQubit has already been called at row, column)
     */
    LogicalQubit *_getLogicalQubit(unsigned int row, unsigned int column);

    /**
     * @brief Get the number of pending gates for the qubits associated with the gate
     *
     * @param gate pointer to the gate
     * @return int the number of pending gates
     */
    int getGateQubitsPendingGatesCount(const Gate *gate);

    /**
     * @brief execute the Rz gates
     *
     * @param rzGates the set of Rz gates to execute
     */
    void _executeRzGates(const std::set<std::pair<int, unsigned int>> &rzGates);

    /**
     * @brief execute the CNOT gates
     *
     * @param cnotGates the set of CNOT gates to execute
     */
    void _executeCNOTGates(const std::set<std::pair<int, unsigned int>> &cnotGates);

    unsigned int _numRows, _numColumns, _numDataQubits, _codeDistance;
    double _physicalQubitErrorRate;
    double (*_rotationErrorModel)(double, double);
    const double _RZ_L_INJECTION_TIME = 2, _RZ_ZZ_INJECTION_TIME = 1, _CLIFFORD_TIME = 2, _CNOT_COMPLETION_TIME = 2, _ROTATION_TIME = 3;
    double _globalTime;
    std::mt19937 _randomNumberGenerator;
    std::vector<Gate *> _gates;
    std::unordered_map<unsigned int, DataQubit *> _dataQubits;
    std::unordered_map<std::pair<unsigned int, unsigned int>, LogicalQubit *, boost::hash<std::pair<unsigned int, unsigned int>>> _logicalQubits;
    std::set<double> _timeSteps;
    std::unordered_map<unsigned int, double> _dataQubitFreeInformationTime;
    std::unordered_map<std::pair<LogicalQubit *, LogicalQubit *>, std::pair<uint8_t, std::vector<LogicalQubit *>>, boost::hash<std::pair<LogicalQubit *, LogicalQubit *>>> _paths;
    std::set<std::pair<int, unsigned int>> _rzGateFrontier, _cnotGateFrontier;
    std::unordered_set<unsigned int> _dataQubitsInFrontier;

    static std::map<std::string, std::vector<double>> _timeTaken;
    static bool _debug;
    static std::ostream *_debugStream;

    static const unsigned int HEATMAP_TIMESTEP;

private:
    /**
     * @brief execute a single Rz gate
     *
     * @param rzGate the Rz gate to execute
     */
    void _executeRzGate(RzGate *rzGate);

    /**
     * @brief process the preparation of the Rz gate
     *
     * @param qubitId the id of the data qubit
     * @param rzGate the Rz gate to prepare
     */
    void _processRzPreparation(unsigned int qubitId, RzGate *rzGate);

    /**
     * @brief execute a single CNOT gate
     *
     * @param cnotGate the CNOT gate to execute
     */
    void _executeCNOTGate(CNOTGate *cnotGate);

    /**
     * @brief update the free information time of the data qubit
     *
     * @param qubitId the id of the data qubit
     * @param time the time relevant to the update
     */
    virtual void _updateDataQubitFreeInformationTime(unsigned int qubitId, double time) = 0;

    /**
     * @brief get the ancilla qubits that can be used for execution (preparation + injection)
     *
     * @param dataQubit the data qubit for which the ancilla qubits are to be found
     * @param isNonClifford whether the gate is non-Clifford
     * @return std::vector<std::pair<AncillaQubit *, AncillaQubit *>> the list of ancilla qubits that can be used for execution (preparation + injection)
     */
    virtual std::vector<std::pair<AncillaQubit *, AncillaQubit *>> _getAncillasOfSingleQubitGate(DataQubit *dataQubit, bool isNonClifford) = 0;

    /**
     * @brief return the required rotation of the qubits (as a bit mask) and the path connecting the two qubits
     *
     * @param cnotGate the CNOT gate for which connection is required
     * @param qubit1 the first qubit
     * @param qubit2 the second qubit
     * @return std::pair<uint8_t, std::vector<LogicalQubit *>> the bitmask of rotation and the path connecting them
     */
    virtual const std::pair<uint8_t, std::vector<LogicalQubit *>> &_connectQubits(CNOTGate *cnotGate, DataQubit *qubit1, DataQubit *qubit2) = 0;

    /**
     * @brief update the paths connecting the qubits (after CNOT execution)
     *
     * @param qubitPair the pair of qubits that were connected
     */
    virtual void _updatePaths(std::pair<LogicalQubit *, LogicalQubit *> qubitPair) = 0;

    /**
     * @brief generate the data qubit indices after compression
     * 
     * @param numColumns the number of columns in the architecture
     * @param numDataQubits the number of data qubits
     * @param compressionFactor the probability of "compression", that is, shifting the qubit to the left or up or down
     * @param seed the seed for the random number generator
     * @return std::tuple<unsigned int, unsigned int, std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>> the number of rows in the compressed architecture, the number of columns in the compressed architecture and the map from data qubit id to the row and column of the qubit
     */
    static std::tuple<unsigned int, unsigned int, std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>> generateDataQubitIndices(unsigned int numColumns, unsigned int numDataQubits, double compressionFactor, double seed = 0);

    std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>> _dataQubitIndices;
};

#endif