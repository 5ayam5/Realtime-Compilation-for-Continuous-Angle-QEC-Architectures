#ifndef __GATE_HPP__
#define __GATE_HPP__

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <boost/functional/hash.hpp>
#include "LogicalQubit.hpp"

/**
 * @brief a generic Quantum Gate class
 *
 */
class Gate
{
public:
    Gate() = delete;

    /**
     * @brief get the name of the gate
     *
     * @return std::string
     */
    std::string getName() const;

    /**
     * @brief Get the start time of the gate
     *
     * @return double
     */
    double getStartTime() const;

    /**
     * @brief Set the start time of the gate
     *
     * @param startTime
     */
    void setStartTime(double startTime);

    /**
     * @brief the less than operator for gates
     *
     * @param other the other gate
     * @return true if this gate has a smaller id
     * @return false otherwise
     */
    bool operator<(const Gate &other) const;

    virtual ~Gate() = default;

protected:
    /**
     * @brief Construct a new Gate object
     *
     * @param name the name of the gate
     */
    Gate(std::string name, unsigned int id);

private:
    std::string _name;
    unsigned int _id;
    double _startTime;
};

/**
 * @brief Rz Gate class
 *
 */
class RzGate : public Gate
{
public:
    /**
     * @brief Construct a new Rz Gate object
     *
     * @param qubit the qubit on which the gate is applied
     * @param theta the rotation angle
     */
    RzGate(unsigned int id, unsigned int qubit, double theta);

    /**
     * @brief get the qubit on which the gate is applied
     *
     * @return unsigned int
     */
    unsigned int getQubit() const;

    /**
     * @brief get the rotation angle
     *
     * @return double
     */
    double getTheta() const;

    /**
     * @brief check if the gate is a non-Clifford gate
     * 
     * @return true if the gate is a non-Clifford gate
     */
    virtual bool isNonClifford() const;

    /**
     * @brief set the rotation angle
     *
     * @param theta the new rotation angle
     */
    virtual void setTheta(double theta);

    /**
     * @brief Get the set of qubits that are preparing the rotation gate for this angle
     *
     * @return std::unordered_set<AncillaQubit *>
     */
    const std::unordered_set<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>, boost::hash<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>>> &getPreparingAncillasCoordinates() const;

    /**
     * @brief add an AncillaQubit to the set of qubits that are preparing the rotation gate for this angle
     *
     * @param ancillaQubit the ancilla qubit to be added
     */
    void addExecutingAncillaQubits(std::pair<AncillaQubit *, AncillaQubit *> ancillaQubits);

    /**
     * @brief remove an AncillaQubit that is being used for preparation of the rotation gate for this angle
     *
     * @param ancillaQubit the AncillaQubit to be removed
     * @param timeOfRemoval the time at which the AncillaQubit is removed
     */
    void removeExecutingAncillaQubits(std::pair<AncillaQubit *, AncillaQubit *> ancillaQubits, double timeOfRemoval);

    /**
     * @brief return true if the gate has succeeded
     *
     * @return true if the gate has succeeded
     * @return false otherwise
     */
    bool hasSucceeded() const;

    /**
     * @brief set the success status of the gate
     *
     * @param succeeded
     */
    void setSucceeded(bool succeeded);

    /**
     * @brief return true if the ancillas are in preparation
     *
     * @return true if ancillas are in preparation
     * @return false otherwise
     */
    bool isInPreparation() const;

    /**
     * @brief Get the expected time it will take to prepare `M_\\theta`
     *
     * @param rotationErrorModel the rotation error model to use
     * @return double
     */
    double getExpectedThetaPrepareTime(unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), std::mt19937 &randomNumberGenerator);

    /**
     * @brief Get the actual time it takes to preare `M_\\theta`
     *
     * @param rotationErrorModel the rotation error model to use
     * @return double
     */
    double getThetaPrepareTime(unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), std::mt19937 &randomNumberGenerator) const;

protected:
    RzGate(unsigned int id, unsigned int qubit, double theta, std::string name);

private:
    unsigned int _qubit;
    double _theta;
    std::unordered_set<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>, boost::hash<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>>> _preparingAncillasCoordinates;
    bool _hasSucceeded;
    bool _inPreparation;

    static std::unordered_map<double, double> _expectedThetaPrepareTimes;
    static const unsigned int NUMBER_OF_SIMULATIONS;
};

/**
 * @brief CNOT Gate class
 *
 */
class CNOTGate : public Gate
{
public:
    /**
     * @brief Construct a new CNOT Gate object
     *
     * @param control the control qubit
     * @param target the target qubit
     */
    CNOTGate(unsigned int id, unsigned int control, unsigned int target);

    /**
     * @brief get the control qubit
     *
     * @return unsigned int
     */
    unsigned int getControl() const;

    /**
     * @brief get the target qubit
     *
     * @return unsigned int
     */
    unsigned int getTarget() const;

private:
    unsigned int _control, _target;
};

class HadamardGate : public RzGate
{
public:
    /**
     * @brief Construct a new Hadamard Gate object
     *
     * @param qubit the qubit on which the gate is applied
     */
    HadamardGate(unsigned int id, unsigned int qubit);

    /**
     * @brief check if the gate is a non-Clifford gate
     * 
     * @return false always since Hadamard gates are Clifford gates
     */
    virtual bool isNonClifford() const override;

    /**
     * @brief set the rotation angle (irrelevant for Hadamard gates)
     */
    virtual void setTheta(double) override;
};

#endif