#include "BaseStarArchitecture.hpp"
#include <cassert>
#include <math.h>
#include <fstream>
#include <chrono>

std::vector<std::pair<unsigned int, unsigned int>> getNeighbourCoordinates(unsigned int row, unsigned int col)
{
    std::vector<std::pair<unsigned int, unsigned int>> neighbours;
    neighbours.push_back({row, col + 1});
    neighbours.push_back({row + 1, col});
    neighbours.push_back({row - 1, col});
    neighbours.push_back({row, col - 1});
    return neighbours;
}

const unsigned int BaseStarArchitecture::HEATMAP_TIMESTEP = 1000;

std::map<std::string, std::vector<double>> BaseStarArchitecture::_timeTaken = std::map<std::string, std::vector<double>>();
bool BaseStarArchitecture::_debug = false;
std::ostream *BaseStarArchitecture::_debugStream = &std::cerr;

BaseStarArchitecture::BaseStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor)
    : _numRows(numRows), _numColumns(numColumns), _numDataQubits(numDataQubits), _codeDistance(codeDistance), _physicalQubitErrorRate(physicalQubitErrorRate), _rotationErrorModel(rotationErrorModel), _globalTime(0), _randomNumberGenerator(seed)
{
    assert(physicalQubitErrorRate >= 0.0f && physicalQubitErrorRate <= 1.0f);

    std::tie(_numRows, _numColumns, _dataQubitIndices) = generateDataQubitIndices(numColumns, numDataQubits, compressionFactor, 2742);
    *_debugStream << "numQubits: " << _numDataQubits << '\n';
    *_debugStream << "numRows: " << _numRows << '\n';
    *_debugStream << "numColumns: " << _numColumns << '\n';
    if (_debug)
    {
        for (auto [qubitId, rowcol] : _dataQubitIndices)
            *_debugStream << "Data qubit " << qubitId << " at row " << rowcol.first << " and column " << rowcol.second << std::endl;
        *_debugStream << std::endl;
    }

    for (unsigned int qubitId = 0; qubitId < numDataQubits; qubitId++)
        _getDataQubit(qubitId);

    _timeSteps.insert(0);
}

unsigned int BaseStarArchitecture::getNumRows() const
{
    return _numRows;
}

unsigned int BaseStarArchitecture::getNumColumns() const
{
    return _numColumns;
}

unsigned int BaseStarArchitecture::getCodeDistance() const
{
    return _codeDistance;
}

const std::vector<Gate *> &BaseStarArchitecture::getGates() const
{
    return _gates;
}

double BaseStarArchitecture::getPhysicalQubitErrorRate() const
{
    return _physicalQubitErrorRate;
}

double BaseStarArchitecture::getGlobalTime() const
{
    return _globalTime;
}

void BaseStarArchitecture::addGates(const std::vector<Gate *> &gates)
{
    for (const auto &gate : gates)
        addGate(gate);
}

void BaseStarArchitecture::addGate(Gate *gate)
{
    if (gate->getName() == "Rz")
    {
        RzGate *rzGate = new RzGate(*static_cast<const RzGate *>(gate));
        _gates.push_back(rzGate);
        unsigned int qubitId = rzGate->getQubit();
        _getDataQubit(qubitId)->pushGate(rzGate);
        if (_debug)
        {
            *_debugStream << "Added Rz gate" << std::endl;
            *_debugStream << "qubit = " << qubitId << std::endl;
            *_debugStream << "theta = " << rzGate->getTheta() << std::endl;
        }
    }
    else if (gate->getName() == "Hadamard")
    {
        HadamardGate *hGate = new HadamardGate(*static_cast<const HadamardGate *>(gate));
        _gates.push_back(hGate);
        unsigned int qubitId = hGate->getQubit();
        _getDataQubit(qubitId)->pushGate(hGate);
        if (_debug)
        {
            *_debugStream << "Added Hadamard gate" << std::endl;
            *_debugStream << "qubit = " << qubitId << std::endl;
        }
    }
    else if (gate->getName() == "CNOT")
    {
        CNOTGate *cnotGate = new CNOTGate(*static_cast<const CNOTGate *>(gate));
        _gates.push_back(cnotGate);
        unsigned int controlId = cnotGate->getControl(), targetId = cnotGate->getTarget();
        _getDataQubit(controlId)->pushGate(cnotGate);
        _getDataQubit(targetId)->pushGate(cnotGate);
        if (_debug)
        {
            *_debugStream << "Added CNOT gate" << std::endl;
            *_debugStream << "control = " << controlId << std::endl;
            *_debugStream << "target = " << targetId << std::endl;
        }
    }
    else
        assert(false);
}

void BaseStarArchitecture::simulate()
{
    while (!step())
        ;
}

void BaseStarArchitecture::logTimeTaken(std::filesystem::path outputFolder)
{
    for (auto [name, times] : _timeTaken)
    {
        std::ofstream file(outputFolder / name);
        double sum = 0;
        *_debugStream << name << ' ';
        if (_debug)
            *_debugStream << std::endl
                          << "    ";
        for (auto time : times)
        {
            if (_debug)
                *_debugStream << time << ' ';
            sum += time;
            file << time << std::endl;
        }
        if (_debug)
            *_debugStream << std::endl
                          << "    ";
        *_debugStream << "average: " << sum / times.size() << '\n';
    }
}

void BaseStarArchitecture::updateQubitsBusyTimes(std::map<unsigned int, std::map<unsigned int, double>> &qubitsBusyTimes, bool ancillas)
{
    for (unsigned int epoch = 0; epoch < _globalTime / HEATMAP_TIMESTEP; epoch++)
        for (unsigned int row = 0; row < _numRows; row++)
            for (unsigned int column = 0; column < _numColumns; column++)
                if (ancillas && _getAncillaQubit(row, column))
                    qubitsBusyTimes[epoch][row * _numColumns + column] += _getAncillaQubit(row, column)->getTotalBusyTime(epoch * HEATMAP_TIMESTEP, (epoch + 1) * HEATMAP_TIMESTEP);
                else if (!ancillas && _getAncillaQubit(row, column) == nullptr)
                {
                    assert(dynamic_cast<DataQubit *>(_getLogicalQubit(row, column)));
                    qubitsBusyTimes[epoch][row * _numColumns + column] += _getLogicalQubit(row, column)->getTotalBusyTime(epoch * HEATMAP_TIMESTEP, (epoch + 1) * HEATMAP_TIMESTEP);
                }
                else
                    qubitsBusyTimes[epoch][row * _numColumns + column] = -1;
}

void BaseStarArchitecture::initStaticStuff(bool debug, std::ostream &debugStream)
{
    _debug = debug;
    _debugStream = &debugStream;
}

void BaseStarArchitecture::resetStaticStuff()
{
    _timeTaken.clear();
}

BaseStarArchitecture::~BaseStarArchitecture()
{
    for (auto logicalQubit : _logicalQubits)
        delete logicalQubit.second;
}

DataQubit *BaseStarArchitecture::_getDataQubit(unsigned int qubitId)
{
    if (_dataQubits.find(qubitId) != _dataQubits.end())
        return _dataQubits[qubitId];

    auto [row, column] = _dataQubitIndices[qubitId];
    DataQubit *dataQubit = new DataQubit(qubitId, row, column);

    _dataQubits.insert({qubitId, dataQubit});
    assert(_logicalQubits.find({row, column}) == _logicalQubits.end());
    _logicalQubits.insert({{row, column}, dataQubit});

    if (_debug)
        *_debugStream << "Initialized data qubit " << qubitId << " at row " << row << " and column " << column << std::endl;
    return dataQubit;
}

AncillaQubit *BaseStarArchitecture::_getAncillaQubit(unsigned int row, unsigned int column)
{
    if (row >= _numRows || column >= _numColumns)
        return nullptr;
    if (_logicalQubits.find({row, column}) == _logicalQubits.end())
    {
        _logicalQubits.insert({{row, column}, new AncillaQubit(row, column)});
        if (_debug)
            *_debugStream << "Initialized ancilla qubit at row " << row << " and column " << column << std::endl;
    }
    return dynamic_cast<AncillaQubit *>(_logicalQubits[{row, column}]);
}

LogicalQubit *BaseStarArchitecture::_getLogicalQubit(unsigned int row, unsigned int column)
{
    if (row >= _numRows || column >= _numColumns)
        return nullptr;
    return _logicalQubits[{row, column}];
}

int BaseStarArchitecture::getGateQubitsPendingGatesCount(const Gate *gate)
{
    auto rzGate = dynamic_cast<const RzGate *>(gate);
    auto cnotGate = dynamic_cast<const CNOTGate *>(gate);
    if (rzGate)
        return _getDataQubit(rzGate->getQubit())->getNumPendingGates();
    else if (cnotGate)
        return std::max(_getDataQubit(cnotGate->getControl())->getNumPendingGates(), _getDataQubit(cnotGate->getTarget())->getNumPendingGates());
    else
    {
        assert(false);
        return 0;
    }
}

void BaseStarArchitecture::_executeRzGates(const std::set<std::pair<int, unsigned int>> &rzGates)
{
    for (auto it = rzGates.begin(); it != rzGates.end();)
    {
        auto [_, qubitId] = *it++;
        auto rzGate = dynamic_cast<RzGate *>(_getDataQubit(qubitId)->topGate());
        _processRzPreparation(qubitId, rzGate);
        if (_debug)
            *_debugStream << "Executing " << rzGate->getName() << " gate on qubit = " << rzGate->getQubit() << std::endl;
        _executeRzGate(rzGate);
    }
}

void BaseStarArchitecture::_executeCNOTGates(const std::set<std::pair<int, unsigned int>> &cnotGates)
{
    for (auto it = cnotGates.begin(); it != cnotGates.end();)
    {
        auto [_, qubitId] = *it++;
        auto cnotGate = static_cast<CNOTGate *>(_getDataQubit(qubitId)->topGate());
        if (_debug)
            *_debugStream << "Attempting to execute CNOT gate between control = " << cnotGate->getControl() << " and target = " << cnotGate->getTarget() << std::endl;
        _executeCNOTGate(cnotGate);
    }
}

void BaseStarArchitecture::_executeRzGate(RzGate *rzGate)
{
    auto clearPreparingAncillas = [this](RzGate *rzGate, std::pair<AncillaQubit *, AncillaQubit *> earliestFreeAncillas = {nullptr, nullptr})
    {
        for (auto it = rzGate->getPreparingAncillasCoordinates().begin(); it != rzGate->getPreparingAncillasCoordinates().end();)
        {
            std::pair<AncillaQubit *, AncillaQubit *> ancillas = {_getAncillaQubit(it->first.first, it->first.second), _getAncillaQubit(it->second.first, it->second.second)};
            ++it;
            if (ancillas != earliestFreeAncillas)
            {
                rzGate->removeExecutingAncillaQubits(ancillas, _globalTime);
                _timeSteps.insert(_globalTime);
                if (_debug)
                    *_debugStream << "    ancilla (" << ancillas.first->getRowNumber() << ", " << ancillas.first->getColumnNumber() << "), (" << (ancillas.second == nullptr ? -1 : ancillas.second->getRowNumber()) << ", " << (ancillas.second == nullptr ? -1 : ancillas.second->getColumnNumber()) << ") freed at " << _globalTime << std::endl;
            }
        }
    };
    rzGate->setStartTime(_globalTime);

    auto preparingAncillasCoordinates = rzGate->getPreparingAncillasCoordinates();
    std::pair<AncillaQubit *, AncillaQubit *> earliestFreeAncillas = {nullptr, nullptr};
    for (auto coordinates : preparingAncillasCoordinates)
    {
        std::pair<AncillaQubit *, AncillaQubit *> ancillas = {_getAncillaQubit(coordinates.first.first, coordinates.first.second), _getAncillaQubit(coordinates.second.first, coordinates.second.second)};
        if (ancillas.first->getBusyUntil() <= _globalTime)
            earliestFreeAncillas = ancillas;
    }
    if (earliestFreeAncillas.first != nullptr)
        clearPreparingAncillas(rzGate, earliestFreeAncillas);

    const auto &dataQubit = _getDataQubit(rzGate->getQubit());
    if (dataQubit->getBusyUntil() > _globalTime)
        return;

    if (rzGate->hasSucceeded())
    {
        _rzGateFrontier.erase({-getGateQubitsPendingGatesCount(rzGate), rzGate->getQubit()});
        _dataQubitsInFrontier.erase(dataQubit->getId());
        clearPreparingAncillas(rzGate);
        dataQubit->removeGate(rzGate);
        _timeSteps.insert(_globalTime);
        _timeTaken[rzGate->getName()].push_back(_globalTime - rzGate->getStartTime());
        if (_debug)
            *_debugStream << "    " << rzGate->getName() << " gate on qubit = " << rzGate->getQubit() << " completed at " << _globalTime << std::endl;
        return;
    }

    if (earliestFreeAncillas.first == nullptr || (earliestFreeAncillas.second != nullptr && (earliestFreeAncillas.second->topGate() != rzGate || earliestFreeAncillas.second->getBusyUntil() > _globalTime)))
        return;
    assert(_globalTime == dataQubit->getBusyUntil() || _globalTime == earliestFreeAncillas.first->getBusyUntil() || earliestFreeAncillas.second != nullptr);

    double injectionEndTime = _globalTime + (earliestFreeAncillas.second == nullptr ? (rzGate->isNonClifford() ? _RZ_ZZ_INJECTION_TIME : _CLIFFORD_TIME) : _RZ_L_INJECTION_TIME);
    dataQubit->setBusy(_globalTime, injectionEndTime);
    earliestFreeAncillas.first->setBusy(_globalTime, injectionEndTime);
    if (earliestFreeAncillas.second != nullptr)
        earliestFreeAncillas.second->setBusy(_globalTime, injectionEndTime);

    rzGate->setTheta(fmod(2 * rzGate->getTheta(), 2 * M_PI));
    rzGate->removeExecutingAncillaQubits(earliestFreeAncillas, injectionEndTime);
    rzGate->setSucceeded(static_cast<double>(_randomNumberGenerator()) / static_cast<double>(_randomNumberGenerator.max()) < 0.5);
    _timeSteps.insert(injectionEndTime);
    _updateDataQubitFreeInformationTime(dataQubit->getId(), _globalTime);
    if (_debug)
        *_debugStream << "    ancilla (" << earliestFreeAncillas.first->getRowNumber() << ", " << earliestFreeAncillas.first->getColumnNumber() << "), (" << (earliestFreeAncillas.second == nullptr ? -1 : earliestFreeAncillas.second->getRowNumber()) << ", " << (earliestFreeAncillas.second == nullptr ? -1 : earliestFreeAncillas.second->getColumnNumber()) << ") injected between " << _globalTime << " and " << injectionEndTime << " (success = " << rzGate->hasSucceeded() << ')' << std::endl;
}

void BaseStarArchitecture::_processRzPreparation(unsigned int qubitId, RzGate *rzGate)
{
    if (!rzGate->isInPreparation() || (rzGate->hasSucceeded() && _getDataQubit(qubitId)->getBusyUntil() == _globalTime))
        return;

    double expectedPrepareTime = rzGate->getExpectedThetaPrepareTime(_codeDistance, _physicalQubitErrorRate, _rotationErrorModel, _randomNumberGenerator);
    double earliestPrepareStartTime = std::max(_dataQubitFreeInformationTime[qubitId], _globalTime - expectedPrepareTime);
    if (earliestPrepareStartTime > _globalTime)
    {
        _timeSteps.insert(earliestPrepareStartTime);
        return;
    }
    if (_debug)
        *_debugStream << rzGate->getName() << " preparation for qubit " << qubitId << std::endl;

    std::set<double> timeSteps;
    double minTime = std::numeric_limits<double>::max();
    auto preparingAncillasCoordinates = rzGate->getPreparingAncillasCoordinates();
    for (auto ancillas : _getAncillasOfSingleQubitGate(_getDataQubit(qubitId), rzGate->isNonClifford()))
    {
        std::pair<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>> coordinatesPair = {{ancillas.first->getRowNumber(), ancillas.first->getColumnNumber()}, {-1, -1}};
        if (ancillas.second)
            coordinatesPair.second = {ancillas.second->getRowNumber(), ancillas.second->getColumnNumber()};
        if (preparingAncillasCoordinates.find(coordinatesPair) != preparingAncillasCoordinates.end())
            continue;
        double actualPrepareTime = rzGate->getThetaPrepareTime(_codeDistance, _physicalQubitErrorRate, _rotationErrorModel, _randomNumberGenerator);
        minTime = std::min(minTime, _globalTime + actualPrepareTime);
        timeSteps.insert(minTime);
        ancillas.first->setBusy(_globalTime, _globalTime + actualPrepareTime, _globalTime + expectedPrepareTime);
        rzGate->addExecutingAncillaQubits(ancillas);

        if (_debug)
            *_debugStream << "    ancilla (" << ancillas.first->getRowNumber() << ", " << ancillas.first->getColumnNumber() << "), (" << (ancillas.second == nullptr ? -1 : ancillas.second->getRowNumber()) << ", " << (ancillas.second == nullptr ? -1 : ancillas.second->getColumnNumber()) << ") prepared between " << _globalTime << " and " << _globalTime + actualPrepareTime << " (expected prepare time was " << expectedPrepareTime << ')' << std::endl;
        if (rzGate->isNonClifford())
            _timeTaken["prepare"].push_back(actualPrepareTime);
    }

    for (auto time : timeSteps)
    {
        _timeSteps.insert(time);
        if (time == minTime)
            break;
    }
}

void BaseStarArchitecture::_executeCNOTGate(CNOTGate *cnotGate)
{
    auto rotateQubit = [this](DataQubit *qubit, AncillaQubit *ancillaQubit) -> bool
    {
        if (ancillaQubit->getBusyUntil() > _globalTime || ancillaQubit->topGate() != qubit->topGate())
            return false;
        qubit->rotate();
        qubit->setBusy(_globalTime, _globalTime + _ROTATION_TIME);
        ancillaQubit->setBusy(_globalTime, _globalTime + _ROTATION_TIME);
        _timeSteps.insert(_globalTime + _ROTATION_TIME);
        if (_debug)
            *_debugStream << "    rotating qubit " << qubit->getId() << " from " << _globalTime << " to " << _globalTime + _ROTATION_TIME << std::endl;
        return true;
    };
    cnotGate->setStartTime(_globalTime);
    DataQubit *control = _getDataQubit(cnotGate->getControl()), *target = _getDataQubit(cnotGate->getTarget());

    auto [rotationMask, path] = _connectQubits(cnotGate, control, target);
    if ((rotationMask >> 1) ^ control->isRotated())
        rotationMask ^= rotateQubit(control, static_cast<AncillaQubit *>(path[1])) << 1;
    if ((rotationMask ^ target->isRotated()) & 1)
        rotationMask ^= rotateQubit(target, static_cast<AncillaQubit *>(path[path.size() - 2]));
    if (rotationMask ^ ((control->isRotated() << 1) | target->isRotated()))
        return;

    double maxTime = 0;
    for (auto qubit : path)
    {
        if (qubit->getBusyUntil() > _globalTime || qubit->topGate() != cnotGate)
            return;
        maxTime = std::max(maxTime, qubit->getBusyUntil());
    }
    double endTime = _globalTime + _CNOT_COMPLETION_TIME;
    _cnotGateFrontier.erase({-getGateQubitsPendingGatesCount(cnotGate), control->getId()});
    _dataQubitsInFrontier.erase(control->getId());
    _dataQubitsInFrontier.erase(target->getId());
    for (auto qubit : path)
    {
        qubit->setBusy(_globalTime, endTime);
        qubit->removeGate(cnotGate);
    }
    _updateDataQubitFreeInformationTime(control->getId(), _globalTime);
    _updateDataQubitFreeInformationTime(target->getId(), _globalTime);
    _timeSteps.insert(endTime);
    _updatePaths(std::minmax(control, target));

    _timeTaken[cnotGate->getName()].push_back(endTime - cnotGate->getStartTime());
    if (_debug)
        *_debugStream << "Executing CNOT gate between control = " << control->getId() << " and target = " << target->getId() << " from " << _globalTime << " to " << endTime << std::endl;
}

std::tuple<unsigned int, unsigned int, std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>>> BaseStarArchitecture::generateDataQubitIndices(unsigned int numColumns, unsigned int numDataQubits, double compressionFactor, double seed)
{
    std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>> dataQubitIndices;
    unsigned int qubitId = 0, rowNum = 1;
    std::mt19937 randomNumberGenerator(seed);
    std::uniform_real_distribution<double> distribution(0, 1);
    unsigned int requiredColumns = 1;

    while (qubitId < numDataQubits)
    {
        // if (distribution(randomNumberGenerator) < compressionFactor / 2)
        //     rowNum ^= 1;
        dataQubitIndices[qubitId++] = {rowNum, 0};

        for (unsigned int columnNum = 2; columnNum < numColumns - 1; columnNum += 2)
        {
            if (qubitId == numDataQubits)
                break;

            double random = distribution(randomNumberGenerator);
            // if (random < compressionFactor / 2)
            //     rowNum ^= 1;
            if (random < compressionFactor)
                columnNum--;
            dataQubitIndices[qubitId++] = {rowNum, columnNum};
            requiredColumns = std::max(requiredColumns, columnNum + 2);
        }
        rowNum |= 1;
        rowNum += 2;
    }

    return {rowNum & ~1, requiredColumns, dataQubitIndices};
}