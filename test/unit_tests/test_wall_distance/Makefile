CXX = g++
CXXFLAGS = -std=c++17 -Wall -I ../../../include -I ~/libraries/include
LDFLAGS = -L ~/libraries/lib
LDLIBS = -lgtest -lgtest_main -pthread

SRC_DIR = ../../../src
EXCLUDE = $(SRC_DIR)/turbobfm.cpp $(SRC_DIR)/turbobfm_grad.cpp
SRC = $(filter-out $(EXCLUDE), $(wildcard $(SRC_DIR)/*.cpp)) main.cpp
OBJ = $(SRC:.cpp=.o)
OBJ_DEBUG = $(SRC:.cpp=.o)

TARGET = main
TARGET_DEBUG = main_debug

all: $(TARGET)

debug: $(TARGET_DEBUG)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

$(TARGET_DEBUG): $(OBJ_DEBUG)
	$(CXX) $(OBJ_DEBUG) $(LDFLAGS) $(LDLIBS) -o $(TARGET_DEBUG)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -O2 -c $< -o $@

%.debug.o: %.cpp
	$(CXX) $(CXXFLAGS) -g -O0 -DDEBUG -c $< -o $@

clean:
	rm -f $(OBJ) $(OBJ_DEBUG) $(TARGET) $(TARGET_DEBUG)