# Pole-zero mismatch simulation. Preferred build: `nix build` (this flake;
# MUSIC tooling from the pinned input) -> result/bin/. This Makefile is the
# interactive fallback, run inside `nix develop`:
#   make -j
# The unmodified MUSIC PulseHistory (and the Constants it reads) are compiled
# from $(MUSIC_TOOLING) into this project's objects; nothing there is edited.

CXX      := g++
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

MUSIC_TOOLING ?= $(HOME)/MUSIC/tooling

INC_DIR   := include
SRC_DIR   := src
MAIN_DIR  := mains
BUILD_DIR := build
BIN_DIR   := bin

CXXFLAGS := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
            $(ROOT_CFLAGS) -I$(INC_DIR) -I$(MUSIC_TOOLING)/include \
            -DMUSIC_HOT_PATH_LOGGING=0 -MMD -MP
LDFLAGS  := $(ROOT_LIBS) -lMinuit -lanalysis-utils -lpthread

SRC_OBJS   := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))
MUSIC_SRCS := PulseHistory Constants
MUSIC_OBJS := $(addprefix $(BUILD_DIR)/music_,$(addsuffix .o,$(MUSIC_SRCS)))
HEADERS    := $(wildcard $(INC_DIR)/*.hpp)

BINS := simulate analyze realwave realdt kernelstudy taustudy figures wavecheck deconvolve optfilter applicability blkernel wavenoise notefigures
BIN_PATHS := $(addprefix $(BIN_DIR)/,$(BINS))

.PHONY: all clean
all: $(BIN_PATHS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/music_%.o: $(MUSIC_TOOLING)/src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR)/main_%.o: $(MAIN_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BIN_DIR)/%: $(BUILD_DIR)/main_%.o $(SRC_OBJS) $(MUSIC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

-include $(wildcard $(BUILD_DIR)/*.d)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
