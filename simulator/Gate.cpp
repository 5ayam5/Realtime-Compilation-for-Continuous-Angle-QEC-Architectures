#include "Gate.hpp"
#include <cassert>
#include <iostream>

template <typename T>
T fast_power(T base, unsigned int exp)
{
    T result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

unsigned int nCr(unsigned int n, unsigned int r)
{
    assert(n >= r);
    if (r > n / 2)
        r = n - r;

    unsigned int result = 1;
    for (unsigned int i = 1; i <= r; i++)
        result *= n - r + i;
    for (unsigned int i = 1; i <= r; i++)
        result /= i;

    return result;
}

Gate::Gate(std::string name, unsigned int id) : _name(name), _id(id), _startTime(-1)
{
}

std::string Gate::getName() const
{
    return _name;
}

double Gate::getStartTime() const
{
    return _startTime;
}

void Gate::setStartTime(double startTime)
{
    if (_startTime == -1)
        _startTime = startTime;
}

bool Gate::operator<(const Gate &other) const
{
    return _id < other._id;
}

std::unordered_map<double, double> RzGate::_expectedThetaPrepareTimes = std::unordered_map<double, double>();
const unsigned int RzGate::NUMBER_OF_SIMULATIONS = 100000;

RzGate::RzGate(unsigned int id, unsigned int qubit, double theta) : RzGate(id, qubit, theta, "Rz")
{
}

unsigned int RzGate::getQubit() const
{
    return _qubit;
}

double RzGate::getTheta() const
{
    return _theta;
}

bool RzGate::isNonClifford() const
{
    return fabs(fmod(_theta, M_PI / 2)) > 1e-6;
}

void RzGate::setTheta(double theta)
{
    _theta = theta;
}

const std::unordered_set<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>, boost::hash<std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>>> &RzGate::getPreparingAncillasCoordinates() const
{
    return _preparingAncillasCoordinates;
}

void RzGate::addExecutingAncillaQubits(std::pair<AncillaQubit *, AncillaQubit *> ancillaQubits)
{
    std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>> coordinatesPair = {{ancillaQubits.first->getRowNumber(), ancillaQubits.first->getColumnNumber()}, {-1, -1}};
    if (ancillaQubits.second != nullptr)
        coordinatesPair.second = {ancillaQubits.second->getRowNumber(), ancillaQubits.second->getColumnNumber()};

    ancillaQubits.first->pushGate(this);
    if (ancillaQubits.second)
        ancillaQubits.second->pushGate(this, ancillaQubits.first);

    assert(_preparingAncillasCoordinates.find(coordinatesPair) == _preparingAncillasCoordinates.end());
    _preparingAncillasCoordinates.insert(coordinatesPair);
}

void RzGate::removeExecutingAncillaQubits(std::pair<AncillaQubit *, AncillaQubit *> ancillaQubits, double timeOfRemoval)
{
    _inPreparation = false;
    std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>> coordinatesPair = {{ancillaQubits.first->getRowNumber(), ancillaQubits.first->getColumnNumber()}, {-1, -1}};
    if (ancillaQubits.second != nullptr)
        coordinatesPair.second = {ancillaQubits.second->getRowNumber(), ancillaQubits.second->getColumnNumber()};

    assert(ancillaQubits.first->topGate() == this);
    ancillaQubits.first->removeGate(this);
    ancillaQubits.first->setFreeAfter(timeOfRemoval);
    if (ancillaQubits.second != nullptr)
        ancillaQubits.second->removeGate(this, ancillaQubits.first);

    assert(_preparingAncillasCoordinates.find(coordinatesPair) != _preparingAncillasCoordinates.end());
    _preparingAncillasCoordinates.erase(coordinatesPair);
}

bool RzGate::hasSucceeded() const
{
    return _hasSucceeded;
}

void RzGate::setSucceeded(bool hasSucceeded)
{
    _inPreparation = true;
    if (!isNonClifford())
        _hasSucceeded = true;
    else
        _hasSucceeded = hasSucceeded;
}

bool RzGate::isInPreparation() const
{
    return _inPreparation;
}

double RzGate::getExpectedThetaPrepareTime(unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), std::mt19937 &randomNumberGenerator)
{
    if (!isNonClifford())
        return 0;
    double probability = (*rotationErrorModel)(physicalQubitErrorRate, _theta);
    if (_expectedThetaPrepareTimes.find(probability) != _expectedThetaPrepareTimes.end())
        return _expectedThetaPrepareTimes[probability];

    double expectedTime = 0;
    for (unsigned int i = 0; i < NUMBER_OF_SIMULATIONS; i++)
        expectedTime += getThetaPrepareTime(codeDistance, physicalQubitErrorRate, rotationErrorModel, randomNumberGenerator);
    expectedTime /= NUMBER_OF_SIMULATIONS;

    return _expectedThetaPrepareTimes[probability] = expectedTime;
}

double RzGate::getThetaPrepareTime(unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), std::mt19937 &randomNumberGenerator) const
{
    double time = 0;
    if (!isNonClifford())
        return time;

    auto circuitSuccessProbability = [physicalQubitErrorRate](unsigned int numberOfOperations)
    {
        double result = 0;
        for (unsigned int i = 0; i < numberOfOperations; i += 2)
            result += nCr(numberOfOperations, i) * fast_power(physicalQubitErrorRate, i) * fast_power(1 - physicalQubitErrorRate, numberOfOperations - i);
        return result;
    };

    /*
        The circuit is as follows:
        GATES --- ROTATION --- ANCILLA PREPARATION --- MEASUREMENT1 --- ANCILLA PREPARATION --- MEASUREMENT2
          4   ---     1    ---         12          ---       4      ---         12          ---       4
        Upto MEASUREMENT1, we can combine all errors and the success probability is stored in `mThetaSuccessUptoMeasurement1` (= `s1`)
        Between MEASUREMENT1 and MEASUREMENT2, we can combine all errors and the success probability is stored in `mThetaSuccessBetweenMeasurement12` (= `s2`)
        The success probability of each group of MEASUREMENTs is stored in `mThetaMeasurementSuccessProbability` (= `sm`)
        The success probability of the whole circuit is stored in `mThetaSuccessProbability`
        `mThetaSuccessProbability` is the probability of getting a `1111` for both measurements under the assumption that each error flips between `1111` and another bit string (and vice versa) with equal probability and the error probability is `physicalQubitErrorRate` for all gates and measurements
    */
    double s1, sm, s2, mThetaSuccessProbability;
    s1 = circuitSuccessProbability(16) * (1 - (*rotationErrorModel)(physicalQubitErrorRate, _theta)) + (1 - circuitSuccessProbability(16)) * (*rotationErrorModel)(physicalQubitErrorRate, _theta);
    sm = fast_power(1 - physicalQubitErrorRate, 4);
    s2 = circuitSuccessProbability(12);
    mThetaSuccessProbability = s1 * sm * (s2 * sm + (1 - s2) * (1 - sm) / 15) + (1 - s1) * (1 - sm) / 15 * ((1 - s2) / 15 * sm + (s2 + 14 * (1 - s2) / 15) * (1 - sm) / 15);
    double logExpansionSuccessProbability = ((codeDistance * codeDistance - 1) / 2) * (log(circuitSuccessProbability(7)) + log(circuitSuccessProbability(5)));

    std::uniform_real_distribution<double> uniformDistribution(0.0f, 1.0f);
    do // expansion success
    {
        do // parallel success of [[4,1,1,2]] preparation
            time += 4.0 / codeDistance;
        while (log(uniformDistribution(randomNumberGenerator)) < fast_power((codeDistance - 1) / 2, 2) * log(1 - mThetaSuccessProbability)); // ((d - 1) / 2)^2 parallel preparations and we loop if all fail

        time += 2.0 / codeDistance;
    } while (logExpansionSuccessProbability < log(uniformDistribution(randomNumberGenerator))); // loop if expansion fails

    return time;
}

RzGate::RzGate(unsigned int id, unsigned int qubit, double theta, std::string name) : Gate(name, id), _qubit(qubit), _theta(theta), _hasSucceeded(false), _inPreparation(true)
{
}

CNOTGate::CNOTGate(unsigned int id, unsigned int control, unsigned int target) : Gate("CNOT", id), _control(control), _target(target)
{
}

unsigned int CNOTGate::getControl() const
{
    return _control;
}

unsigned int CNOTGate::getTarget() const
{
    return _target;
}

HadamardGate::HadamardGate(unsigned int id, unsigned int qubit) : RzGate(id, qubit, INFINITY, "Hadamard")
{
}

bool HadamardGate::isNonClifford() const
{
    return false;
}

void HadamardGate::setTheta(double)
{
}