#ifndef __LOGICALQUBIT_HPP__
#define __LOGICALQUBIT_HPP__

#include <list>
#include <vector>

class Gate;
class AncillaQubit;

class LogicalQubit
{
public:
    LogicalQubit() = delete;

    virtual ~LogicalQubit() = default;

    /**
     * @brief get the row number of the logical qubit
     *
     * @return unsigned int
     */
    unsigned int getRowNumber() const;

    /**
     * @brief get the column number of the logical qubit
     *
     * @return unsigned int
     */
    unsigned int getColumnNumber() const;

    /**
     * @brief get the coordinates for qubits along the Z edges of the logical qubit
     *
     * @return std::vector<std::pair<unsigned int, unsigned int>> the neighbouring coordinates
     */
    std::vector<std::pair<unsigned int, unsigned int>> getZNeighbours() const;

    /**
     * @brief get the coordinates for qubits along the X edges of the logical qubit
     *
     * @return std::vector<std::pair<unsigned int, unsigned int>> the neighbouring coordinates
     */
    std::vector<std::pair<unsigned int, unsigned int>> getXNeighbours() const;

    /**
     * @brief check if the qubit is currently busy
     *
     * @param time the time at which we need to check
     * @return true if the qubit is busy
     * @return false if the qubit is free
     */
    bool isBusy(double time) const;

    /**
     * @brief get the latest time until which the qubit is busy
     *
     * @return double
     */
    double getBusyUntil() const;

    /**
     * @brief set the qubit to be busy until a certain time
     *
     * @param startTime the time from which the qubit is busy
     * @param endTime the time until which the qubit is busy
     * @param expectedEndTime the time until which the qubit is expected to be busy
     */
    void setBusy(double startTime, double endTime, double expectedEndTime = -1);

    /**
     * @brief Get the total time the qubit is busy between two times
     *
     * @param startTime the start time of the interval
     * @param endTime the end time of the interval
     * @return double the total busy time
     */
    double getTotalBusyTime(double startTime, double endTime) const;

    /**
     * @brief Get the time until which the qubit is expected to be busy
     *
     * @param queryTime the time at which the query is made
     *
     * @return double
     */
    virtual double getExpectedBusyUntil(double queryTime) = 0;

    /**
     * @brief push a gate to the queue of operations
     *
     * @param gate the gate to push
     * @param associatedAncilla the ancilla qubit associated with the gate
     */
    void pushGate(Gate *gate, AncillaQubit *associatedAncilla = nullptr);

    /**
     * @brief peek the next gate to execute from the queue
     *
     * @return constant pointer to the gate and nullptr if the queue is empty
     */
    Gate *topGate() const;

    /**
     * @brief pop the next gate from the queue
     *
     * @param gate the gate to remove
     * @param associatedAncilla the ancilla qubit associated with the gate
     */
    void removeGate(Gate *gate, AncillaQubit *associatedAncilla = nullptr);

    /**
     * @brief get the number of pending gates
     *
     * @return unsigned int the number of pending gates
     */
    unsigned int getNumPendingGates() const;

protected:
    LogicalQubit(unsigned int rowNumber, unsigned int columnNumber);

    unsigned int _rowNumber, _columnNumber;
    bool _isRotated;
    std::vector<std::pair<double, double>> _busyTimes;
    double _expectedBusyUntil;
    std::list<std::pair<Gate *, AncillaQubit *>> _gateOperations;
};

class AncillaQubit : public LogicalQubit
{
public:
    AncillaQubit(unsigned int rowNumber, unsigned int columnNumber);

    /**
     * @brief set the ancilla free after a certain time (if it is no longer needed)
     *
     * @param time the time after which the ancilla is free
     */
    void setFreeAfter(double time);

    /**
     * @brief Get the time until which the qubit is expected to be busy
     *
     * @param queryTime the time at which the query is made
     *
     * @return double
     */
    double getExpectedBusyUntil(double queryTime) override;
};

/**
 * @brief Data Qubit class
 *
 */
class DataQubit : public LogicalQubit
{
public:
    /**
     * @brief Construct a new Logical Qubit object
     *
     * @param id the id of the logical qubit
     * @param rowNumber the row number of the logical qubit
     * @param columnNumber the column number of the logical qubit
     */
    DataQubit(unsigned int id, unsigned int rowNumber, unsigned int columnNumber);

    /**
     * @brief get the id of the logical qubit
     *
     * @return unsigned int
     */
    unsigned int getId() const;

    /**
     * @brief rotate the logical qubit
     *
     */
    void rotate();

    /**
     * @brief check if the logical qubit is rotated
     *
     * @return true if the logical qubit is rotated
     * @return false otherwise
     */
    bool isRotated() const;

    /**
     * @brief Get the time until which the qubit is expected to be busy
     *
     * @return double
     */
    double getExpectedBusyUntil(double) override;

private:
    unsigned int _id;
};

#endif