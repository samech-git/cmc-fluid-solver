CUDA_INSTALL_PATH = /usr/local/cuda
NETCDF_INSTALL_PATH = /usr

CXX := g++
CC := gcc
LINK := g++ -fPIC
NVCC := nvcc -ccbin /usr/bin -arch=compute_10 -code=sm_10

# Includes
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include -I$(NETCDF_INSTALL_PATH)/include -I../FluidSolver2D

# Common flags
COMMONFLAGS += $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)

# Libs
LIB_CUDA :=  -L$(CUDA_INSTALL_PATH)/lib64 -lcudart
LIB_NETCDF := -L$(NETCDF_INSTALL_PATH)/lib64 -lnetcdf

OBJS = Grid2D.cpp.o AdiSolver3D.cpp.o AdiSolver3D.cu.o Grid3D.cpp.o \
       TimeLayer3D.cu.o Solver3D.cpp.o FluidSolver3D.cpp.o LinuxIO.cpp.o
TARGET = ../../bin/Release/FluidSolver3D
LINKLINE = $(LINK) -o $(TARGET) $(OBJS) $(LIB_CUDA) $(LIB_NETCDF) -lrt

.SUFFIXES: .cpp .cu .o

%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

Grid2D.cpp.o: ../FluidSolver2D/Grid2D.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

LinuxIO.cpp.o: ../Common/LinuxIO.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: $(OBJS) Makefile
	$(LINKLINE)

clean:
	rm $(TARGET) *.o
