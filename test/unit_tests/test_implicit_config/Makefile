# Compiler
CXX = g++
CXXFLAGS = -std=c++20 -Wall -I ../../../include/ -I ~/libraries/include -Xclang -fopenmp -I /opt/homebrew/opt/libomp/include

# Libraries
LDFLAGS = -L ~/libraries/lib -L /opt/homebrew/opt/libomp/lib -lomp
LDLIBS = -lgtest -lgtest_main -pthread

# Sources and objects
SRC = main.cpp ../../../src/config.cpp
OBJ = $(SRC:.cpp=.o)

# Target
TARGET = main

# Default rule
all: $(TARGET)

# Link step
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET)
