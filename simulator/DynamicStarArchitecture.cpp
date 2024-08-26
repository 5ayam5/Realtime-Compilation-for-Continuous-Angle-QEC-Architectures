#include "DynamicStarArchitecture.hpp"
#include "DSU.hpp"

DynamicStarArchitecture::DynamicStarArchitecture(unsigned int numRows, unsigned int numColumns, unsigned int numDataQubits, unsigned int codeDistance, double physicalQubitErrorRate, double (*rotationErrorModel)(double, double), unsigned int seed, double mstComputationFrequency, double compressionFactor)
    : BaseStarArchitecture(numRows, numColumns, numDataQubits, codeDistance, physicalQubitErrorRate, rotationErrorModel, seed, compressionFactor), _mstComputationFrequency(mstComputationFrequency)
{
    computeMST(-_MST_COMPUTATION_TIME);
}

bool DynamicStarArchitecture::step()
{
    if (_rzGateFrontier.empty() && _cnotGateFrontier.empty())
    {
        bool done = true;
        for (auto qubit : _dataQubits)
            if (qubit.second->topGate() != nullptr || qubit.second->getBusyUntil() > _globalTime)
                done = false;
        if (done)
            return true;
    }

    _globalTime = *_timeSteps.begin();
    _timeSteps.erase(_timeSteps.begin());
    if (_debug)
        *_debugStream << "Executing DynamicStarArchitecture::step() at time = " << _globalTime << std::endl;

    if (fmod(_globalTime, _mstComputationFrequency) == 0)
        computeMST(_globalTime);

    for (auto [qubitId, dataQubit] : _dataQubits)
    {
        if (_dataQubitsInFrontier.find(qubitId) != _dataQubitsInFrontier.end())
            continue;
        auto const gate = dataQubit->topGate();
        if (!gate)
            continue;
        if (gate->getName() == "Rz" || gate->getName() == "Hadamard")
        {
            if (dataQubit->getBusyUntil() > _globalTime)
                continue;
            _rzGateFrontier.insert({-getGateQubitsPendingGatesCount(gate), qubitId});
            _dataQubitsInFrontier.insert(qubitId);
        }
        else if (gate->getName() == "CNOT")
        {
            auto cnotGate = static_cast<CNOTGate *>(gate);
            DataQubit *control = _getDataQubit(cnotGate->getControl()), *target = _getDataQubit(cnotGate->getTarget());
            if (control->getBusyUntil() > _globalTime || target->getBusyUntil() > _globalTime)
                continue;
            if (control->topGate() != gate || target->topGate() != gate)
                continue;
            _cnotGateFrontier.insert({-getGateQubitsPendingGatesCount(gate), control->getId()});
            _dataQubitsInFrontier.insert(control->getId());
            _dataQubitsInFrontier.insert(target->getId());
        }
        else
        {
            if (_debug)
                *_debugStream << "Gate name " << gate->getName() << " not supported" << std::endl;
            assert(false);
        }
    }

    _executeCNOTGates(_cnotGateFrontier);
    _executeRzGates(_rzGateFrontier);

    if (_debug)
        *_debugStream << std::endl;

    return false;
}

void DynamicStarArchitecture::_updateDataQubitFreeInformationTime(unsigned int qubitId, double time)
{
    _dataQubitFreeInformationTime[qubitId] = time;
}

std::vector<std::pair<AncillaQubit *, AncillaQubit *>> DynamicStarArchitecture::_getAncillasOfSingleQubitGate(DataQubit *dataQubit, bool isNonClifford)
{
    std::vector<std::pair<AncillaQubit *, AncillaQubit *>> ancillas;

    auto neighbours = dataQubit->getZNeighbours();
    for (auto neighbour : neighbours)
    {
        auto ancillaQubit = _getAncillaQubit(neighbour.first, neighbour.second);
        if (ancillaQubit && ancillaQubit->topGate() == nullptr && ancillaQubit->getBusyUntil() <= _globalTime)
            ancillas.push_back({ancillaQubit, nullptr});
    }

    neighbours = dataQubit->getXNeighbours();
    for (auto neighbour : neighbours)
    {
        auto ancillaQubit = _getAncillaQubit(neighbour.first, neighbour.second);
        if (ancillaQubit)
        {
            if (!isNonClifford)
            {
                if (ancillaQubit->topGate() == nullptr && ancillaQubit->getBusyUntil() <= _globalTime)
                    ancillas.push_back({ancillaQubit, nullptr});
            }
            else
            {
                auto nextNeighbours = dataQubit->isRotated() ? ancillaQubit->getXNeighbours() : ancillaQubit->getZNeighbours();
                for (auto nextNeighbour : nextNeighbours)
                {
                    auto nextAncillaQubit = _getAncillaQubit(nextNeighbour.first, nextNeighbour.second);
                    if (nextAncillaQubit && nextAncillaQubit->topGate() == nullptr && nextAncillaQubit->getBusyUntil() <= _globalTime)
                        ancillas.push_back({nextAncillaQubit, ancillaQubit});
                }
            }
        }
    }

    return ancillas;
}

std::pair<double, std::vector<LogicalQubit *>> DynamicStarArchitecture::_getPath(std::unordered_set<std::pair<AncillaQubit *, AncillaQubit *>, boost::hash<std::pair<AncillaQubit *, AncillaQubit *>>> tree, DataQubit *qubit1, AncillaQubit *neighbour1, DataQubit *qubit2, AncillaQubit *neighbour2)
{
    std::vector<LogicalQubit *> path;
    path.push_back(qubit1);
    std::unordered_set<AncillaQubit *> visited;
    std::function<double(AncillaQubit *)> dfs;
    dfs = [this, &visited, &path, &tree, &dfs, neighbour2](AncillaQubit *qubit) -> double
    {
        visited.insert(qubit);
        path.push_back(qubit);
        if (qubit == neighbour2)
            return qubit->getBusyUntil();

        for (auto neighbour : getNeighbourCoordinates(qubit->getRowNumber(), qubit->getColumnNumber()))
        {
            AncillaQubit *neighbourQubit = _getAncillaQubit(neighbour.first, neighbour.second);
            if (neighbourQubit != nullptr && tree.find(std::minmax(qubit, neighbourQubit)) != tree.end() && visited.find(neighbourQubit) == visited.end())
            {
                double ret = dfs(neighbourQubit);
                if (ret != -1)
                    return std::max(ret, qubit->getBusyUntil());
            }
        }

        path.pop_back();
        return -1;
    };
    auto ret = dfs(neighbour1);
    if (ret == -1)
        return {-1, {}};
    path.push_back(qubit2);
    return {ret, path};
}

const std::pair<uint8_t, std::vector<LogicalQubit *>> &DynamicStarArchitecture::_connectQubits(CNOTGate *cnotGate, DataQubit *qubit1, DataQubit *qubit2)
{
    if (_paths.find(std::minmax(qubit1, qubit2)) != _paths.end())
        return _paths[std::minmax(qubit1, qubit2)];
    _msts.erase(_msts.begin(), --_msts.upper_bound(_globalTime));

    double earliestStartTime = std::max(_globalTime, std::max(qubit1->getBusyUntil(), qubit2->getBusyUntil()));
    std::map<double, std::pair<uint8_t, std::vector<LogicalQubit *>>> paths;
    for (bool controlRotated : {false, true})
        for (bool targetRotated : {false, true})
        {
            std::vector<std::pair<unsigned int, unsigned int>> controlNeighbours, targetNeighbours;
            if (controlRotated)
                controlNeighbours = qubit1->getXNeighbours();
            else
                controlNeighbours = qubit1->getZNeighbours();
            if (targetRotated)
                targetNeighbours = qubit2->getZNeighbours();
            else
                targetNeighbours = qubit2->getXNeighbours();

            for (auto controlNeighbour : controlNeighbours)
                for (auto targetNeighbour : targetNeighbours)
                {
                    auto controlAncillaQubit = _getAncillaQubit(controlNeighbour.first, controlNeighbour.second), targetAncillaQubit = _getAncillaQubit(targetNeighbour.first, targetNeighbour.second);
                    if (controlAncillaQubit != nullptr && targetAncillaQubit != nullptr)
                    {
                        auto [startTime, path] = _getPath(_msts.begin()->second, qubit1, controlAncillaQubit, qubit2, targetAncillaQubit);
                        if (startTime == -1)
                            continue;
                        startTime = std::max(startTime, earliestStartTime);
                        if (controlRotated)
                        {
                            double rotationStartTime = std::max(_globalTime, std::max(qubit1->getBusyUntil(), controlAncillaQubit->getExpectedBusyUntil(_globalTime)));
                            startTime = std::max(startTime, rotationStartTime + _ROTATION_TIME);
                        }
                        if (targetRotated)
                        {
                            double rotationStartTime = std::max(_globalTime, std::max(qubit2->getBusyUntil(), targetAncillaQubit->getExpectedBusyUntil(_globalTime)));
                            startTime = std::max(startTime, rotationStartTime + _ROTATION_TIME);
                        }
                        if (paths.find(startTime) == paths.end())
                            paths[startTime] = {(controlRotated ^ qubit1->isRotated()) << 1 | (targetRotated ^ qubit2->isRotated()), path};
                    }
                }
        }

    assert(!paths.empty());
    auto &path = _paths[std::minmax(qubit1, qubit2)] = paths.begin()->second;
    if (_debug)
        *_debugStream << "    Routing CNOT between " << qubit1->getId() << " and " << qubit2->getId() << std::endl;
    for (auto qubit : path.second)
    {
        AncillaQubit *ancillaQubit = dynamic_cast<AncillaQubit *>(qubit);
        if (ancillaQubit)
        {
            if (_debug)
                *_debugStream << "        ancilla at (" << ancillaQubit->getRowNumber() << ", " << ancillaQubit->getColumnNumber() << ") reserved" << std::endl;
            ancillaQubit->pushGate(cnotGate);
        }
        else
            assert(static_cast<DataQubit *>(qubit) == qubit1 || static_cast<DataQubit *>(qubit) == qubit2);
    }
    return _paths[std::minmax(qubit1, qubit2)] = path;
}

void DynamicStarArchitecture::_updatePaths(std::pair<LogicalQubit *, LogicalQubit *> qubitPair)
{
    assert(_paths.find(qubitPair) != _paths.end());
    _paths.erase(qubitPair);
}

void DynamicStarArchitecture::computeMST(double startTime)
{
    auto processQubit = [this](DSU<AncillaQubit *> &dsu, std::set<std::pair<double, std::pair<AncillaQubit *, AncillaQubit *>>> &edges, AncillaQubit *ancillaQubit)
    {
        dsu.make_set(ancillaQubit);
        double busyTime = ancillaQubit->getTotalBusyTime(_globalTime - _MST_COMPUTATION_TIME, _globalTime);
        auto addEdge = [this, &edges, ancillaQubit, busyTime](std::pair<unsigned int, unsigned int> neighbour)
        {
            AncillaQubit *neighbourQubit = _getAncillaQubit(neighbour.first, neighbour.second);
            if (neighbourQubit)
                edges.insert({std::max(busyTime, neighbourQubit->getTotalBusyTime(_globalTime - _MST_COMPUTATION_TIME, _globalTime)), std::minmax(ancillaQubit, neighbourQubit)});
        };
        for (auto neighbour : ancillaQubit->getXNeighbours())
            addEdge(neighbour);
        for (auto neighbour : ancillaQubit->getZNeighbours())
            addEdge(neighbour);
    };

    auto time = std::chrono::high_resolution_clock::now();
    DSU<AncillaQubit *> dsu;
    std::set<std::pair<double, std::pair<AncillaQubit *, AncillaQubit *>>> edges;

    for (unsigned int row = 0; row < _numRows; row++)
        for (unsigned int column = 0; column < _numColumns; column++)
        {
            auto ancillaQubit = _getAncillaQubit(row, column);
            if (ancillaQubit)
                processQubit(dsu, edges, ancillaQubit);
        }

    std::unordered_set<std::pair<AncillaQubit *, AncillaQubit *>, boost::hash<std::pair<AncillaQubit *, AncillaQubit *>>> mst;
    for (auto edge : edges)
    {
        auto qubit1 = edge.second.first, qubit2 = edge.second.second;
        if (dsu.find_set(qubit1) != dsu.find_set(qubit2))
        {
            dsu.union_set(qubit1, qubit2);
            mst.insert(std::minmax(qubit1, qubit2));
        }
    }
    _timeTaken["mst"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - time).count());

    if (_debug)
        *_debugStream << "    MST computed between time = " << startTime << " and " << startTime + _MST_COMPUTATION_TIME << std::endl;

    _msts[startTime + _MST_COMPUTATION_TIME] = mst;
    _timeSteps.insert(startTime + _mstComputationFrequency);
}