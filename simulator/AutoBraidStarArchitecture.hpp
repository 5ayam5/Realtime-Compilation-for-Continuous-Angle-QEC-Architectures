#ifndef __AUTOBRAID_STAR_ARCHITECTURE_HPP__
#define __AUTOBRAID_STAR_ARCHITECTURE_HPP__

#include "StaticStarArchitecture.hpp"

class AutoBraidStarArchitecture : public StaticStarArchitecture
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
    AutoBraidStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor = 0);

    /**
     * @brief simulate the architecture until the next gate
     *
     * @return true if there are no more gates to execute
     */
    bool step() override;

private:
    /**
     * @brief return the required rotation of the qubits (as a bit mask) and the path connecting the two qubits
     *
     * @param qubit1 the first qubit
     * @param qubit2 the second qubit
     * @return std::pair<uint8_t, std::vector<LogicalQubit *>> the bitmask of rotation and the path connecting them
     */
    virtual const std::pair<uint8_t, std::vector<LogicalQubit *>> &_connectQubits(CNOTGate *, DataQubit *qubit1, DataQubit *qubit2) override;

    /**
     * @brief construct the paths for as many CNOTs as possible using AutoBraid (https://dl.acm.org/doi/10.1145/3466752.3480072) and insert them into the frontier
     * 
     * @param cnotGates the list of CNOT gates in the current layer
     */
    void _constructPaths(std::vector<CNOTGate *> &cnotGates);
};

#endif