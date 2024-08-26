#include "StaticStarArchitecture.hpp"
#include <queue>

StaticStarArchitecture::StaticStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor)
    : BaseStarArchitecture(numRows, numColumns, numDataQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, compressionFactor)
{
}

bool StaticStarArchitecture::step()
{
    _globalTime = *_timeSteps.begin();
    _timeSteps.erase(_timeSteps.begin());
    if (_debug)
        *_debugStream << "Executing StaticStarArchitecture::step() at time = " << _globalTime << std::endl;

    if (!_rzGateFrontier.empty() || !_cnotGateFrontier.empty())
    {
        _executeRzGates(_rzGateFrontier);
        _executeCNOTGates(_cnotGateFrontier);
        return false;
    }
    if (!_timeSteps.empty())
        return false;
    else
    {
        bool done = true;
        for (auto qubit : _dataQubits)
            if (qubit.second->topGate() != nullptr)
                done = false;
        if (done)
            return true;
    }

    _paths.clear();
    _freeQubits.clear();
    for (unsigned int row = 0; row < _numRows; row++)
        for (unsigned int column = 0; column < _numColumns; column++)
        {
            auto ancillaQubit = _getAncillaQubit(row, column);
            if (ancillaQubit)
                _freeQubits.insert(ancillaQubit);
        }

    for (auto [qubitId, dataQubit] : _dataQubits)
    {
        auto const gate = dataQubit->topGate();
        auto rzGate = dynamic_cast<RzGate *>(gate);
        auto cnotGate = dynamic_cast<CNOTGate *>(gate);
        if (rzGate)
        {
            assert(cnotGate == nullptr);
            auto dataQubit = _getDataQubit(rzGate->getQubit());
            if (_dataQubitsInFrontier.find(dataQubit->getId()) != _dataQubitsInFrontier.end())
                continue;
            assert(dataQubit->topGate() == rzGate);

            auto ancillas = _getAncillasOfSingleQubitGate(dataQubit, rzGate->isNonClifford());
            if (_freeQubits.find(ancillas[0].first) == _freeQubits.end() || (ancillas[0].second && _freeQubits.find(ancillas[0].second) == _freeQubits.end()))
                continue;

            if (_debug)
                *_debugStream << "Adding " << rzGate->getName() << " gate (" << rzGate << ") to frontier for qubit " << rzGate->getQubit() << std::endl;

            _dataQubitFreeInformationTime[dataQubit->getId()] = _globalTime;
            _rzGateFrontier.insert({-getGateQubitsPendingGatesCount(rzGate), dataQubit->getId()});
            _dataQubitsInFrontier.insert(rzGate->getQubit());
            _freeQubits.erase(ancillas[0].first);
            _freeQubits.erase(ancillas[0].second);
        }
        else if (cnotGate)
        {
            assert(cnotGate->getControl() < _numDataQubits);
            assert(cnotGate->getTarget() < _numDataQubits);
            auto control = _getDataQubit(cnotGate->getControl()), target = _getDataQubit(cnotGate->getTarget());
            if (_dataQubitsInFrontier.find(control->getId()) != _dataQubitsInFrontier.end() || _dataQubitsInFrontier.find(target->getId()) != _dataQubitsInFrontier.end() || control->topGate() != cnotGate || target->topGate() != cnotGate)
                continue;

            auto [rotationMask, path] = _connectQubits(cnotGate, control, target);
            if (rotationMask == (uint8_t)-1)
                continue;
            for (auto qubit : path)
            {
                _freeQubits.erase(qubit);
                if (dynamic_cast<AncillaQubit *>(qubit) != nullptr)
                    qubit->pushGate(cnotGate);
            }

            if (_debug)
                *_debugStream << "Adding CNOT gate to frontier for control " << cnotGate->getControl() << " and target " << cnotGate->getTarget() << std::endl;

            _cnotGateFrontier.insert({-getGateQubitsPendingGatesCount(cnotGate), control->getId()});
            _dataQubitsInFrontier.insert(control->getId());
            _dataQubitsInFrontier.insert(target->getId());
        }
        else
            assert(false);
    }

    _executeRzGates(_rzGateFrontier);
    _executeCNOTGates(_cnotGateFrontier);

    return false;
}

std::vector<std::pair<AncillaQubit *, AncillaQubit *>> StaticStarArchitecture::_getAncillasOfSingleQubitGate(DataQubit *dataQubit, bool isNonClifford)
{
    if (isNonClifford)
        return std::vector<std::pair<AncillaQubit *, AncillaQubit *>>({{_getAncillaQubit(dataQubit->getRowNumber() ^ 1, dataQubit->getColumnNumber() + 1), _getAncillaQubit(dataQubit->getRowNumber() ^ 1, dataQubit->getColumnNumber())}});
    else
        return std::vector<std::pair<AncillaQubit *, AncillaQubit *>>({{_getAncillaQubit(dataQubit->getRowNumber() ^ 1, dataQubit->getColumnNumber()), nullptr}});
}

std::vector<LogicalQubit *> StaticStarArchitecture::_findShortestPath(DataQubit *qubit1, DataQubit *qubit2)
{
    std::vector<LogicalQubit *> path;
    std::unordered_set<LogicalQubit *> visited;
    std::queue<std::pair<std::vector<LogicalQubit *>, LogicalQubit *>> q;
    visited.insert(qubit1);
    for (auto neighbour : getNeighbourCoordinates(qubit1->getRowNumber(), qubit1->getColumnNumber()))
    {
        AncillaQubit *neighbourQubit = _getAncillaQubit(neighbour.first, neighbour.second);
        if (neighbourQubit)
        {
            q.push({{qubit1, neighbourQubit}, neighbourQubit});
            visited.insert(neighbourQubit);
        }
    }
    while (!q.empty())
    {
        auto [currentPath, currentQubit] = q.front();
        q.pop();
        if (currentQubit == qubit2)
        {
            path = currentPath;
            break;
        }
        for (auto neighbour : getNeighbourCoordinates(currentQubit->getRowNumber(), currentQubit->getColumnNumber()))
        {
            LogicalQubit *neighbourQubit = _getLogicalQubit(neighbour.first, neighbour.second);
            if ((neighbourQubit == qubit2 || _freeQubits.find(neighbourQubit) != _freeQubits.end()) && visited.find(neighbourQubit) == visited.end())
            {
                visited.insert(neighbourQubit);
                auto newPath = currentPath;
                newPath.push_back(neighbourQubit);
                q.push({newPath, neighbourQubit});
            }
        }
    }
    return path;
}

void StaticStarArchitecture::_updateDataQubitFreeInformationTime(unsigned int qubitId, double)
{
    _dataQubitFreeInformationTime[qubitId] = _getDataQubit(qubitId)->getBusyUntil();
}

const std::pair<uint8_t, std::vector<LogicalQubit *>> &StaticStarArchitecture::_connectQubits(CNOTGate *, DataQubit *qubit1, DataQubit *qubit2)
{
    if (_paths.find(std::minmax(qubit1, qubit2)) != _paths.end())
        return _paths[std::minmax(qubit1, qubit2)];


    auto path = _findShortestPath(qubit1, qubit2);
    if (path.empty())
        return _paths[std::minmax(qubit1, qubit2)] = {(uint8_t)-1, {}};

    bool controlRotated = false, targetRotated = false;
    auto neighbours = qubit1->getXNeighbours();
    for (auto neighbour : neighbours)
        if (_getAncillaQubit(neighbour.first, neighbour.second) == path[1])
            controlRotated = true;
    neighbours = qubit2->getZNeighbours();
    for (auto neighbour : neighbours)
        if (_getAncillaQubit(neighbour.first, neighbour.second) == path[path.size() - 2])
            targetRotated = true;

    if (_debug)
    {
        *_debugStream << "Path between " << qubit1->getId() << " and " << qubit2->getId() << ":" << std::endl;
        for (auto qubit : path)
            *_debugStream << '(' << qubit->getRowNumber() << ", " << qubit->getColumnNumber() << ") ";
        *_debugStream << std::endl;
    }

    return _paths[std::minmax(qubit1, qubit2)] = {((controlRotated ^ qubit1->isRotated()) << 1) | (targetRotated ^ qubit2->isRotated()), path};
}

void StaticStarArchitecture::_updatePaths(std::pair<LogicalQubit *, LogicalQubit *> qubitPair)
{
    auto rotateQubit = [this](DataQubit *qubit, AncillaQubit *ancillaQubit)
    {
        double startTime = std::max(qubit->getBusyUntil(), ancillaQubit->getBusyUntil());
        qubit->rotate();
        qubit->setBusy(startTime, startTime + _ROTATION_TIME);
        ancillaQubit->setBusy(startTime, startTime + _ROTATION_TIME);
        _timeSteps.insert(startTime + _ROTATION_TIME);
        if (_debug)
            *_debugStream << "    rotating qubit " << qubit->getId() << " from " << startTime << " to " << startTime + _ROTATION_TIME << std::endl;
        return true;
    };
    assert(_paths.find(qubitPair) != _paths.end());
    auto [rotationMask, path] = _paths[qubitPair];

    auto control = static_cast<DataQubit *>(path[0]), target = static_cast<DataQubit *>(path[path.size() - 1]);
    if ((rotationMask >> 1) ^ control->isRotated())
        rotateQubit(control, static_cast<AncillaQubit *>(path[1]));
    if ((rotationMask ^ target->isRotated()) & 1)
        rotateQubit(target, static_cast<AncillaQubit *>(path[path.size() - 2]));

    _paths.erase(qubitPair);
}