#include "AutoBraidStarArchitecture.hpp"
#include <stack>
#include <utility>

AutoBraidStarArchitecture::AutoBraidStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double compressionFactor)
    : StaticStarArchitecture(numRows, numColumns, numDataQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, compressionFactor)
{
}

bool AutoBraidStarArchitecture::step()
{
    _globalTime = *_timeSteps.begin();
    _timeSteps.erase(_timeSteps.begin());
    if (_debug)
        *_debugStream << "Executing AutoBraidArchitecture::step() at time = " << _globalTime << std::endl;

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
        auto rzGate = dynamic_cast<RzGate *>(dataQubit->topGate());
        if (rzGate)
        {
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
    }

    std::vector<CNOTGate *> cnotGates;
    for (auto [qubitId, dataQubit] : _dataQubits)
    {
        auto cnotGate = dynamic_cast<CNOTGate *>(dataQubit->topGate());
        if (cnotGate)
        {
            assert(cnotGate->getControl() < _numDataQubits);
            assert(cnotGate->getTarget() < _numDataQubits);
            auto control = _getDataQubit(cnotGate->getControl()), target = _getDataQubit(cnotGate->getTarget());
            if (control != dataQubit || _dataQubitsInFrontier.find(control->getId()) != _dataQubitsInFrontier.end() || _dataQubitsInFrontier.find(target->getId()) != _dataQubitsInFrontier.end() || control->topGate() != cnotGate || target->topGate() != cnotGate)
                continue;
            cnotGates.push_back(cnotGate);
        }
    }
    _constructPaths(cnotGates);

    _executeRzGates(_rzGateFrontier);
    _executeCNOTGates(_cnotGateFrontier);

    return false;
}

const std::pair<uint8_t, std::vector<LogicalQubit *>> &AutoBraidStarArchitecture::_connectQubits(CNOTGate *, DataQubit *qubit1, DataQubit *qubit2)
{
    assert(_paths.find(std::minmax(qubit1, qubit2)) != _paths.end());
    return _paths.at(std::minmax(qubit1, qubit2));
}

void AutoBraidStarArchitecture::_constructPaths(std::vector<CNOTGate *> &cnotGates)
{
    auto getBoudingBox = [this](CNOTGate *cnotGate)
    {
        auto control = _getDataQubit(cnotGate->getControl()), target = _getDataQubit(cnotGate->getTarget());
        return std::make_pair(std::minmax(control->getRowNumber(), target->getRowNumber()), std::minmax(control->getColumnNumber(), target->getColumnNumber()));
    };

    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> adjList;
    for (uint32_t i = 0; i < cnotGates.size(); i++)
        for (uint32_t j = i + 1; j < cnotGates.size(); j++)
        {
            auto cnot1 = cnotGates[i], cnot2 = cnotGates[j];
            auto boundingBox1 = getBoudingBox(cnot1), boundingBox2 = getBoudingBox(cnot2);
            if (boundingBox1.first.first <= boundingBox2.first.second && boundingBox1.first.second >= boundingBox2.first.first && boundingBox1.second.first <= boundingBox2.second.second && boundingBox1.second.second >= boundingBox2.second.first)
            {
                adjList[i].insert(j);
                adjList[j].insert(i);
            }
        }

    std::unordered_set<uint32_t> cnotGatesSet;
    for (uint32_t i = 0; i < cnotGates.size(); i++)
        cnotGatesSet.insert(i);

    auto getArea = [getBoudingBox](CNOTGate *cnotGate) {
        auto boundingBox = getBoudingBox(cnotGate);
        return (boundingBox.first.second - boundingBox.first.first) * (boundingBox.second.second - boundingBox.second.first);
    };

    std::stack<uint32_t> stack;
    while (!cnotGatesSet.empty())
    {
        uint32_t bestIdx = *cnotGatesSet.begin(), bestArea = getArea(cnotGates[bestIdx]);
        for (auto idx : cnotGatesSet)
        {
            auto area = getArea(cnotGates[idx]);
            if (area > bestArea || (area == bestArea && adjList[idx].size() > adjList[bestIdx].size()))
            {
                bestArea = area;
                bestIdx = idx;
            }
        }
        stack.push(bestIdx);
        cnotGatesSet.erase(bestIdx);
        for (auto neighbour : adjList[bestIdx])
            adjList.at(neighbour).erase(bestIdx);
        adjList.erase(bestIdx);
    }

    while (!stack.empty())
    {
        auto idx = stack.top();
        stack.pop();
        auto cnotGate = cnotGates[idx];
        auto control = _getDataQubit(cnotGate->getControl()), target = _getDataQubit(cnotGate->getTarget());

        assert(_paths.find(std::minmax(control, target)) == _paths.end());
        auto path = _findShortestPath(control, target);
        if (path.empty())
            continue;

        bool controlRotated = false, targetRotated = false;
        auto neighbours = control->getXNeighbours();
        for (auto neighbour : neighbours)
            if (_getAncillaQubit(neighbour.first, neighbour.second) == path[1])
                controlRotated = true;
        neighbours = target->getZNeighbours();
        for (auto neighbour : neighbours)
            if (_getAncillaQubit(neighbour.first, neighbour.second) == path[path.size() - 2])
                targetRotated = true;
        _paths[std::minmax(control, target)] = {((controlRotated ^ control->isRotated()) << 1) | (targetRotated ^ target->isRotated()), path};

        if (_debug)
        {
            *_debugStream << "Path between " << control->getId() << " and " << target->getId() << ":" << std::endl;
            for (auto qubit : path)
                *_debugStream << '(' << qubit->getRowNumber() << ", " << qubit->getColumnNumber() << ") ";
            *_debugStream << std::endl;
        }

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
}