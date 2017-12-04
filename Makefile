SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

CC = g++
CYTHONC = cython3
PYTHONINCD = -I/usr/include/python3.5m/
# TODO environment variables
INCD = -I$(HOME)/networkit/include/
LIBD = -L$(HOME)/networkit

LIB = -lNetworKit
SRCDIR = ./src
OBJDIR = ./.build
SRC := $(addprefix $(SRCDIR)/,BioMaxentStress.cpp, BioMaxentStress.cpp, IDGPOptimizerOld.cpp, BioMaxentStressOldOld.cpp, _MOBi.cpp)
# OBJ := $(addprefix $(OBJDIR)/,SRC:.cpp=.o)
OBJ := $(addprefix $(OBJDIR)/,BioMaxentStress.o, BioMaxentStressOld.o, IDGPOptimizerOld.o, BioMaxentStressOldOld.o, _MOBi.o)

CFLAGS = -std=c++11 -Wall -Werror -fopenmp -fPIC -O3 $(INCD)
LFLAGS = -Wl,--whole-archive $(HOME)/networkit/libNetworKit.a -Wl,--no-whole-archive
CYTHONFLAGS = --cplus -Werror -3 -I$(HOME)/networkit/networkit/

all: $(OBJ) ./lib/libMOBi.a ./lib/_MOBi.so

# TODO clean this up
#./lib/_MOBi.so: $(SRCDIR)/_MOBi.cpp ./lib/libMOBi.a
./lib/_MOBi.so: $(SRCDIR)/_MOBi.cpp ./src/BioMaxentStressOldOld.cpp ./src/BioMaxentStressOldOld.h ./src/IDGPOptimizerOld.cpp ./src/IDGPOptimizerOld.h
	# $(CC) -shared -pthread -fwrapv -fno-strict-aliasing $(CFLAGS) $(LFLAGS) $(PYTHONINCD) -o 
	#$(CC) -shared -pthread -fwrapv -fno-strict-aliasing -std=c++11 -Wall -fopenmp -fPIC -O3 $(INCD) $(PYTHONINCD) -o ./lib/_MOBi.so ./src/_MOBi.cpp ./src/BioMaxentStress.cpp  $(LFLAGS) -Wno-sign-compare -Werror
	$(CC) -shared -pthread -fwrapv -fno-strict-aliasing -std=c++11 -Wall -fopenmp -fPIC -O3 -ggdb $(INCD) $(PYTHONINCD) -o ./lib/_MOBi.so ./src/_MOBi.cpp ./src/BioMaxentStressOldOld.cpp  $(LFLAGS) -Wno-sign-compare

# TODO make cython optional i.e. ship the cpp source
$(SRCDIR)/_MOBi.cpp: $(SRCDIR)/_MOBi.pyx
	$(CYTHONC) $(CYTHONFLAGS) -o $(SRCDIR)/_MOBi.cpp $(SRCDIR)/_MOBi.pyx

./lib/libMOBi.a: ./src/BioMaxentStress.cpp ./src/BioMaxentStressOld.cpp ./src/IDGPOptimizerOld.cpp
	ar rcs -o $@ $^

$(OBJDIR)%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	$(RM) $(OBJDIR)/*.o $(SRCDIR)/_MOBi.cpp ./lib/libMOBi.a ./lib/_MOBi.so
