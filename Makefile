# ==== Forked from Jason Lau's Makefile ====
# File:             Makefile
# Description:      Makefile for the problem specified by file path
# Author:           Jason Lau <i@dotkrnl.com>
# Date:             2015-03-15
# License:          Apache
# Feel free to contact me if there's any questions
#

MAIN ?= MeshSimplification

ifeq ($(OS), Windows_NT)
	ARGV ?= 
else
	ARGV ?= 
endif

SRCDIR ?= src
OBJDIR ?= bin

CXX ?= g++
CXXFLAGS ?= -O2 -c -Wall -std=c++11
LDFLAGS ?=
E ?= @echo 

CXXFILES := $(wildcard $(SRCDIR)/*.cpp)
EXE_CXXFILES := $(wildcard $(SRCDIR)/main_*.cpp)
LIB_CXXFILES := $(wildcard $(SRCDIR)/lib_*.cpp)
LIB_OBJFILES := $(LIB_CXXFILES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HFILES := $(wildcard $(SRCDIR)/*.h)

EXECUTABLES := $(EXE_CXXFILES:$(SRCDIR)/main_%.cpp=%)

ifeq ($(OS), Windows_NT)
	MAIN := $(MAIN).exe
	EXECUTABLES := $(EXECUTABLES:%=%.exe)
 	CXXFLAGS += -DWIN32
endif

all: $(EXECUTABLES)
debug: CXXFLAGS += -DDEBUG -g
debug: $(EXECUTABLES)

run: $(MAIN)
	./$(MAIN) $(ARGV)

clean:
	$(E) + Removing files
ifeq ($(OS), Windows_NT)
	del $(EXECUTABLES) $(OBJDIR)\*
else
	rm $(EXECUTABLES) $(OBJDIR)/*
endif

.PRECIOUS: $(OBJDIR)/%.o
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HFILES)
	$(E) + Compiling $<
	$(CXX) -o $@ -c $< $(CXXFLAGS)

ifeq ($(OS), Windows_NT)
%.exe: $(OBJDIR)/main_%.o $(LIB_OBJFILES)
else
%: $(OBJDIR)/main_%.o $(LIB_OBJFILES)
endif
	$(E) + Linking $@
	$(CXX) -o $@ $< $(LIB_OBJFILES) $(LDFLAGS)

