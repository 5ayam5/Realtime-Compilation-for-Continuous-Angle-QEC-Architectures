#ifndef __STATICSTARARCHITECTURE_HPP__
#define __STATICSTARARCHITECTURE_HPP__

#include "BaseStarArchitecture.hpp"
#include <unordered_set>

class StaticStarArchitecture : public BaseStarArchitecture
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
    StaticStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor = 0);

    /**
     * @brief simulate the architecture until the next gate
     *
     * @return true if there are no more gates to execute
     */
    virtual bool step() override;

protected:
    /**
     * @brief get the ancilla qubits that can be used for preparation
     *
     * @param dataQubit the data qubit for which the ancilla qubits are to be found
     * @param isNonClifford whether the gate is non-Clifford
     * @return std::vector<AncillaQubit *> the list of ancilla qubits that can be used for preparation
     */
    std::vector<std::pair<AncillaQubit *, AncillaQubit *>> _getAncillasOfSingleQubitGate(DataQubit *dataQubit, bool isNonClifford) override;

    /**
     * @brief BFS to find the shortest path between two qubits
     * 
     */
    std::vector<LogicalQubit *> _findShortestPath(DataQubit *qubit1, DataQubit *qubit2);

    std::unordered_set<LogicalQubit *> _freeQubits;

private:
    /**
     * @brief update the free information time of the data qubit
     *
     * @param qubitId the id of the data qubit
     * @param time the time relevant to the update
     */
    void _updateDataQubitFreeInformationTime(unsigned int qubitId, double time) override;

    /**
     * @brief return the required rotation of the qubits (as a bit mask) and the path connecting the two qubits
     *
     * @param qubit1 the first qubit
     * @param qubit2 the second qubit
     * @return std::pair<uint8_t, std::vector<LogicalQubit *>> the bitmask of rotation and the path connecting them
     */
    virtual const std::pair<uint8_t, std::vector<LogicalQubit *>> &_connectQubits(CNOTGate *, DataQubit *qubit1, DataQubit *qubit2) override;

    /**
     * @brief update the paths connecting the qubits (after CNOT execution)
     *
     * @param qubitPair the pair of qubits that were connected
     */
    void _updatePaths(std::pair<LogicalQubit *, LogicalQubit *> qubitPair) override;
};

#endif