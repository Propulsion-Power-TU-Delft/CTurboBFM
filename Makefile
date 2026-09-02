# ==========================================================
# Directories
# ==========================================================
SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

PROFILE_LIB_PATH := /opt/homebrew/opt/gperftools/lib
PROFILE_INCLUDE_PATH := /opt/homebrew/opt/gperftools/include

OMP_LIB_PATH := /opt/homebrew/opt/libomp/lib
OMP_INCLUDE_PATH := /opt/homebrew/opt/libomp/include

# ==========================================================
# Compiler and flags
# ==========================================================
CXX := g++
CXXFLAGS := -std=c++20 -Wall -I$(INC_DIR) -I$(SRC_DIR) -Xclang -fopenmp -I$(OMP_INCLUDE_PATH)

LDFLAGS  := -L$(OMP_LIB_PATH) -lomp

# Build-type specific flags
DEBUG_FLAGS   := -g -O0
RELEASE_FLAGS := -O3 -DNDEBUG -march=native
PROFILE_FLAGS := -g -O3 -DNDEBUG -march=native -I$(PROFILE_INCLUDE_PATH)

# Default build type
BUILD_TYPE := debug

# ==========================================================
# Sources and objects
# ==========================================================
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)

# Define main and grad entry points
MAIN_SRC := $(SRC_DIR)/turbobfm.cpp
GRAD_SRC := $(SRC_DIR)/turbobfm_grad.cpp

# Common (shared) sources
COMMON_SOURCES := $(filter-out $(MAIN_SRC) $(GRAD_SRC), $(SOURCES))

# Object lists
COMMON_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(COMMON_SOURCES))
MAIN_OBJECTS   := $(COMMON_OBJECTS) $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(MAIN_SRC))
GRAD_OBJECTS   := $(COMMON_OBJECTS) $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(GRAD_SRC))

# Output binaries
MAIN_TARGET := $(BIN_DIR)/turbobfm
GRAD_TARGET := $(BIN_DIR)/turbobfm_grad

# ==========================================================
# Default rules
# ==========================================================
all: $(MAIN_TARGET)

debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: all

release: CXXFLAGS += $(RELEASE_FLAGS)
release: all

profile: CXXFLAGS += $(PROFILE_FLAGS)
profile: LDFLAGS += -L$(PROFILE_LIB_PATH) -lprofiler
profile: all

# ==========================================================
# Linking rules
# ==========================================================
$(MAIN_TARGET): $(MAIN_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(MAIN_OBJECTS) $(LDFLAGS) -o $(MAIN_TARGET)
	@echo "✅ Built $(MAIN_TARGET)"

$(GRAD_TARGET): $(GRAD_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(GRAD_OBJECTS) $(LDFLAGS) -o $(GRAD_TARGET)
	@echo "✅ Built $(GRAD_TARGET)"

# ==========================================================
# Compilation rule
# ==========================================================
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==========================================================
# Utility rules
# ==========================================================
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "🧹 Cleaned build directories."

# Optional shortcut targets
main: $(MAIN_TARGET)
grad: $(GRAD_TARGET)

.PHONY: all clean debug release profile main grad