#ifndef __DYNAMICSTARARCHITECTURE_HPP__
#define __DYNAMICSTARARCHITECTURE_HPP__

#include "BaseStarArchitecture.hpp"
#include <unordered_set>

class DynamicStarArchitecture : public BaseStarArchitecture
{
public:
    /**
     * @brief Construct a new Patch Star Architecture object
     *
     * @param numRows the number of rows in the architecture
     * @param numColumns the number of columns in the architecture
     * @param codeDistance the surface code distance
     * @param physicalQubitErrorRate the physical qubit error rate
     * @param rotationErrorModel the rotation error model
     * @param seed the seed for the random number generator
     * @param patchNumberRows the number of rows in each patch
     * @param patchNumberColumns the number of columns in each patch
     * @param compressionFactor the probability of "compression", that is, shifting the qubit to the left or up or down
     */
    DynamicStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double mstComputationFrequency, double compressionFactor = 0);

    /**
     * @brief simulate the architecture until the next gate
     *
     * @return true if there are no more gates to execute
     */
    bool step() override;

private:
    /**
     * @brief update the free information time of the data qubit
     *
     * @param qubitId the id of the data qubit
     * @param time the time relevant to the update
     */
    void _updateDataQubitFreeInformationTime(unsigned int qubitId, double time) override;

    /**
     * @brief get the ancilla qubits that can be used for preparation
     *
     * @param dataQubit the data qubit for which the ancilla qubits are to be found
     * @param isNonClifford whether the gate is non-Clifford
     * @return std::vector<AncillaQubit *> the list of ancilla qubits that can be used for preparation
     */
    std::vector<std::pair<AncillaQubit *, AncillaQubit *>> _getAncillasOfSingleQubitGate(DataQubit *dataQubit, bool isNonClifford) override;

    /**
     * @brief get the path connecting the two qubits
     *
     * @param tree the tree of the grid from which the path is to be searched
     * @param qubit1 the first qubit
     * @param qubit2 the second qubit
     * @return std::pair<double, std::vector<LogicalQubit *>> the earliest time that they can be connected and the path connecting them
     */
    std::pair<double, std::vector<LogicalQubit *>> _getPath(std::unordered_set<std::pair<AncillaQubit *, AncillaQubit *>, boost::hash<std::pair<AncillaQubit *, AncillaQubit *>>> tree, DataQubit *qubit1, AncillaQubit *neighbour1, DataQubit *qubit2, AncillaQubit *neighbour2);

    /**
     * @brief return the required rotation of the qubits (as a bit mask) and the path connecting the two qubits
     *
     * @param cnotGate the CNOT gate for which connection is required
     * @param qubit1 the first qubit
     * @param qubit2 the second qubit
     * @return std::pair<uint8_t, std::vector<LogicalQubit *>> the bitmask of rotation and the path connecting them
     */
    const std::pair<uint8_t, std::vector<LogicalQubit *>> &_connectQubits(CNOTGate *cnotGate, DataQubit *qubit1, DataQubit *qubit2) override;

    /**
     * @brief update the paths connecting the qubits (after CNOT execution)
     *
     * @param qubitPair the pair of qubits that were connected
     */
    void _updatePaths(std::pair<LogicalQubit *, LogicalQubit *> qubitPair) override;

    void computeMST(double startTime);

    std::map<double, std::unordered_set<std::pair<AncillaQubit *, AncillaQubit *>, boost::hash<std::pair<AncillaQubit *, AncillaQubit *>>>> _msts;
    const double _MST_COMPUTATION_TIME = 100;
    double _mstComputationFrequency;
};

#endif