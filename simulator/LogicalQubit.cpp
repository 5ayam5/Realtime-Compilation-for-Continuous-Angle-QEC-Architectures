#include "LogicalQubit.hpp"
#include "Gate.hpp"
#include <algorithm>

LogicalQubit::LogicalQubit(unsigned int rowNumber, unsigned int columnNumber)
    : _rowNumber(rowNumber), _columnNumber(columnNumber), _isRotated(false), _busyTimes({{0, 0}}), _expectedBusyUntil(0)
{
}

unsigned int LogicalQubit::getRowNumber() const
{
    return _rowNumber;
}

unsigned int LogicalQubit::getColumnNumber() const
{
    return _columnNumber;
}

std::vector<std::pair<unsigned int, unsigned int>> LogicalQubit::getZNeighbours() const
{
    if (_isRotated)
        return std::vector<std::pair<unsigned int, unsigned int>>({{_rowNumber - 1, _columnNumber}, {_rowNumber + 1, _columnNumber}});
    else
        return std::vector<std::pair<unsigned int, unsigned int>>({{_rowNumber, _columnNumber - 1}, {_rowNumber, _columnNumber + 1}});
}

std::vector<std::pair<unsigned int, unsigned int>> LogicalQubit::getXNeighbours() const
{
    if (_isRotated)
        return std::vector<std::pair<unsigned int, unsigned int>>({{_rowNumber, _columnNumber - 1}, {_rowNumber, _columnNumber + 1}});
    else
        return std::vector<std::pair<unsigned int, unsigned int>>({{_rowNumber - 1, _columnNumber}, {_rowNumber + 1, _columnNumber}});
}

bool LogicalQubit::isBusy(double time) const
{
    auto it = std::lower_bound(_busyTimes.begin(), _busyTimes.end(), time, [](const std::pair<double, double> &busyTime, double time)
                               { return busyTime.first < time; });
    if (it == _busyTimes.end())
        return false;
    return time < it->second;
}

double LogicalQubit::getBusyUntil() const
{
    return _busyTimes.back().second;
}

void LogicalQubit::setBusy(double startTime, double endTime, double expectedEndTime)
{
    _busyTimes.push_back({startTime, endTime});
    if (expectedEndTime == -1)
        _expectedBusyUntil = endTime;
    else
        _expectedBusyUntil = expectedEndTime;
}

double LogicalQubit::getTotalBusyTime(double startTime, double endTime) const
{
    double totalBusyTime = 0;
    for (auto it = _busyTimes.begin(); it != _busyTimes.end(); ++it)
    {
        if (it->second <= startTime)
            continue;
        if (it->first >= endTime)
            break;
        totalBusyTime += std::min(it->second, endTime) - std::max(it->first, startTime);
    }
    return totalBusyTime;
}

void LogicalQubit::pushGate(Gate *gate, AncillaQubit *associatedAncilla)
{
    _gateOperations.push_back({gate, associatedAncilla});
}

Gate *LogicalQubit::topGate() const
{
    if (_gateOperations.empty())
        return nullptr;
    return _gateOperations.front().first;
}

void LogicalQubit::removeGate(Gate *gate, AncillaQubit *associatedAncilla)
{
    auto it = std::find(_gateOperations.begin(), _gateOperations.end(), std::make_pair(gate, associatedAncilla));
    assert(it != _gateOperations.end());
    _gateOperations.erase(it);
}

unsigned int LogicalQubit::getNumPendingGates() const
{
    return _gateOperations.size();
}

AncillaQubit::AncillaQubit(unsigned int rowNumber, unsigned int columnNumber)
    : LogicalQubit(rowNumber, columnNumber)
{
}

void AncillaQubit::setFreeAfter(double time)
{
    if (_busyTimes.back().second <= time)
        return;
    else if (_busyTimes.back().first > time)
        _busyTimes.pop_back();
    else
        _busyTimes.back().second = time;
    _expectedBusyUntil = std::min(time, _busyTimes.back().second);
}

double AncillaQubit::getExpectedBusyUntil(double queryTime)
{
    if (getBusyUntil() <= queryTime)
        return _expectedBusyUntil = getBusyUntil();
    return _expectedBusyUntil;
}

DataQubit::DataQubit(unsigned int id, unsigned int rowNumber, unsigned int columnNumber)
    : LogicalQubit(rowNumber, columnNumber), _id(id)
{
}

unsigned int DataQubit::getId() const
{
    return _id;
}

void DataQubit::rotate()
{
    _isRotated = !_isRotated;
}

bool DataQubit::isRotated() const
{
    return _isRotated;
}

double DataQubit::getExpectedBusyUntil(double)
{
    return getBusyUntil();
}