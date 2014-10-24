CUDA_INSTALL_PATH ?= /usr/local/cuda

CXX := g++ -g
CC := gcc
LINK := g++ -fPIC
NVCC  := $(CUDA_INSTALL_PATH)/bin/nvcc -ccbin /usr/bin

# Includes
INCLUDES = -I. -I$(CUDA_INSTALL_PATH)/include

# Common flags
CXXFLAGS = -Wall -std=c++0x -O2
NVCCFLAGS = -Xcompiler -O2 -gencode=arch=compute_35,code=sm_35
DEFINES = 
LIBS = -lgomp -fopenmp -lfftw3 -lm -llapack -lblas
COMMONFLAGS = $(LIBS) $(DEFINES) $(INCLUDES)
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)

LIB_CUDA := -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcufft -lcublas -lcufftw

# Project Specific things
OBJS = main.cpp deconvolver.o imager/image.o timer.o
TARGET = deconvolve
LINKLINE = $(LINK) $(DEFINES) $(INCLUDES) $(LIB_CUDA) -o $(TARGET) $(OBJS) $(LIBS)


.SUFFIXES: .c .cpp .cu .o

%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(LINKLINE)


run: $(TARGET)
	./$(TARGET) imager/vis.out clean.out imager/img.out

clean:
	-rm *.o
	-rm clean.out
	-rm clean.out.png
	-rm $(TARGET)