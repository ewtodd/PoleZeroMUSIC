# Figures for doc/kernel_note.tex. Preferred build: `nix build` -> result/bin/.
# This Makefile is the interactive fallback, run inside `nix develop`: make -j

CXX         := g++
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

INC_DIR   := include
SRC_DIR   := src
MAIN_DIR  := mains
BUILD_DIR := build
BIN_DIR   := bin

CXXFLAGS := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
            $(ROOT_CFLAGS) -I$(INC_DIR) -MMD -MP
LDFLAGS  := $(ROOT_LIBS) -lanalysis-utils

HEADERS  := $(wildcard $(INC_DIR)/*.hpp)
SRC_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))
BINS     := figures
BIN_PATHS := $(addprefix $(BIN_DIR)/,$(BINS))

.PHONY: all clean
all: $(BIN_PATHS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/main_%.o: $(MAIN_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN_DIR)/%: $(BUILD_DIR)/main_%.o $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

-include $(wildcard $(BUILD_DIR)/*.d)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
