CXX_FLAGS=-std=c++17 -Wall -Wextra -pedantic -Ofast -g -DNDEBUG -I/opt/homebrew/opt/boost/include # replace with path to Boost on your device
CXX=g++

all: sim

SIM_SOURCES := $(wildcard simulator/*.cpp)
SIM_OBJECTS := $(SIM_SOURCES:.cpp=.o)
SIM_HEADERS := $(wildcard simulator/*.hpp)

$(SIM_OBJECTS): %.o: %.cpp $(SIM_HEADERS)
	$(CXX) $(CXX_FLAGS) -c $< -o $@

sim: $(SIM_OBJECTS) $(SIM_HEADERS)
	$(CXX) $(CXX_FLAGS) $(SIM_OBJECTS) -o $@

clean:
	rm -f **/*.o sim

.PHONY: clean
