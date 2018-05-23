SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o

CC = g++
CYTHONC = cython
PYTHONINCD = -I/usr/include/python3.6m/
# TODO environment variables
# TODO this is horrible awefulness
INCD = -I$(HOME)/Projects/networkit/networkit/
# LIBD = -L$(HOME)/Projects/networkit
LIBD = /usr/lib/python3.6/site-packages

LIB = -lNetworKit
SRCDIR = ./src
OBJDIR = ./.build
# SRC := $(addprefix $(SRCDIR)/,BioMaxentStress.cpp, BioMaxentStress.cpp, IDGPOptimizerOld.cpp, BioMaxentStressOldOld.cpp, _MOBi.cpp)
SRC := $(addprefix $(SRCDIR)/,BioMaxentStress.cpp, BioMaxentStress.cpp, BioMaxentStress.cpp, _MOBi.cpp)
# OBJ := $(addprefix $(OBJDIR)/,SRC:.cpp=.o)
OBJ := $(addprefix $(OBJDIR)/,BioMaxentStress.o, _MOBi.o)

# TODO mkdir for lib
CFLAGS = -std=c++11 -Wall -Werror -fopenmp -fPIC -O3 $(INCD)
# LFLAGS = -Wl,--whole-archive $(HOME)/networkit/libNetworKit.a -Wl,--no-whole-archive
LFLAGS = -Wl,--whole-archive $(LIBD)/_NetworKit.cpython-36m-x86_64-linux-gnu.so -Wl,--no-whole-archive
# CYTHONFLAGS = --cplus -Werror -3 -I$(HOME)/networkit/networkit/
CYTHONFLAGS = --cplus -Werror -3 -I$(LIBD)/networkit/networkit/

# TODO create output folders

# all: $(OBJ) ./lib/libMOBi.a ./lib/_MOBi.so
all: $(OBJ) ./lib/_MOBi.so

# TODO clean this up
#./lib/_MOBi.so: $(SRCDIR)/_MOBi.cpp ./lib/libMOBi.a
./lib/_MOBi.so: $(SRCDIR)/_MOBi.cpp ./src/BioMaxentStress.cpp ./src/BioMaxentStress.h
	# $(CC) -shared -pthread -fwrapv -fno-strict-aliasing $(CFLAGS) $(LFLAGS) $(PYTHONINCD) -o 
	#$(CC) -shared -pthread -fwrapv -fno-strict-aliasing -std=c++11 -Wall -fopenmp -fPIC -O3 $(INCD) $(PYTHONINCD) -o ./lib/_MOBi.so ./src/_MOBi.cpp ./src/BioMaxentStress.cpp  $(LFLAGS) -Wno-sign-compare -Werror
	$(CC) -shared -pthread -fwrapv -fno-strict-aliasing -std=c++11 -Wall -fopenmp -fPIC -O3 -ggdb $(INCD) $(PYTHONINCD) -o ./lib/_MOBi.so ./src/_MOBi.cpp ./src/BioMaxentStress.cpp  $(LFLAGS) -Wno-sign-compare

# TODO make cython optional i.e. ship the cpp source
$(SRCDIR)/_MOBi.cpp: $(SRCDIR)/_MOBi.pyx
	$(CYTHONC) $(CYTHONFLAGS) -o $(SRCDIR)/_MOBi.cpp $(SRCDIR)/_MOBi.pyx

# ./lib/libMOBi.a: ./src/BioMaxentStress.cpp ./src/BioMaxentStressOld.cpp ./src/IDGPOptimizerOld.cpp
# 	ar rcs -o $@ $^

$(OBJDIR)%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

.PHONY : clean
clean:
	$(RM) $(OBJDIR)/*.o $(SRCDIR)/_MOBi.cpp ./lib/libMOBi.a ./lib/_MOBi.so
