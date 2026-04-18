CXX := g++

TARGET := fluid_sim
BINDIR := bin
OBJDIR := build/obj

SRC := \
	main.cpp \
	aux.cpp \
	core-sim-functions.cpp \
	display-functions.cpp \
	initializations.cpp \
	imgui/imgui.cpp \
	imgui/imgui_draw.cpp \
	imgui/imgui_tables.cpp \
	imgui/imgui_widgets.cpp \
	imgui-sfml/imgui-SFML.cpp

OBJ := $(SRC:%.cpp=$(OBJDIR)/%.o)
DEP := $(OBJ:.o=.d)

SFML_CFLAGS ?= $(shell pkg-config --cflags sfml-graphics 2>/dev/null)
SFML_LIBS ?= $(shell pkg-config --libs sfml-graphics 2>/dev/null)

ifeq ($(strip $(SFML_LIBS)),)
	SFML_LIBS := -lsfml-graphics -lsfml-window -lsfml-system
endif

CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra
CPPFLAGS += -I. -Iimgui -Iimgui-sfml -DIMGUI_USER_CONFIG='"imconfig-SFML.h"' $(SFML_CFLAGS)
LDFLAGS += $(SFML_LIBS) -lGL

.PHONY: all build run clean help

all: build

build: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: build
	./$(BINDIR)/$(TARGET)

clean:
	rm -rf build $(BINDIR)

help:
	@echo "Targets:"
	@echo "  make build   - Compile the project"
	@echo "  make run     - Compile and run the project"
	@echo "  make clean   - Remove build artifacts"

-include $(DEP)
