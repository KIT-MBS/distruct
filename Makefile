SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

CC = clang
INCD = -I$(HOME)/networkit/include/
LIBD = -L$(HOME)/networkit
LIB = -lNetworKit

SRCDIR = ./src
OBJDIR = ./.build
SRC := $(addprefix $(SRCDIR)/,BioMaxentStress.cpp)
# OBJ := $(addprefix $(OBJDIR)/,SRC:.cpp=.o)
OBJ := $(addprefix $(OBJDIR)/,BioMaxentStress.o)

CFLAGS = -std=c++11 -Wall -fopenmp -O3 $(INCD)

all: $(OBJ)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $<

$(OBJ): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	rm .build/*.o
