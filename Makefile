CUDA     := /usr/local/cuda

SM50     := -gencode=arch=compute_50,code=\"sm_50,compute_50\"
TARGET   := $(SM50)
NVCC     := $(CUDA)/bin/nvcc
LINK     := g++ -m64 -fPIC 
NVFLAGS  := -m64 --ptxas-options -v
INCLUDES := -I. -I$(CUDA)/include -I/usr/local/cuda/include
LIB      := -L$(CUDA)/lib64 -L/usr/local/cuda/lib -L/usr/local/cuda/lib64 -lcuda -lcudart -lntl -lgmp -lm

all: gpuMul.o Base.o DeviceManager.o NTT.o kernel.o 
	$(LINK)  -o all $^ $(LIB)
gpuMul.o:gpuMul.cu
	$(NVCC) -G -g -c $(NVFLAGS) $(TARGET) $(INCLUDES) $(LIB) gpuMul.cu
Base.o:Base.cu Base.h
	$(NVCC) -G -g -c $(NVFLAGS) $(TARGET) $(INCLUDES) $(LIB) Base.cu
DeviceManager.o:DeviceManager.cu DeviceManager.h
	$(NVCC) -G -g -c $(NVFLAGS) $(TARGET) $(INCLUDES) $(LIB) DeviceManager.cu
NTT.o:NTT.cu NTT.h
	$(NVCC) -G -g -c $(NVFLAGS) $(TARGET) $(INCLUDES) $(LIB) NTT.cu
kernel.o:kernel.cu kernel.h
	$(NVCC) -G -g -c $(NVFLAGS) $(TARGET) $(INCLUDES) $(LIB) kernel.cu
clean:
	rm -f main *.o *.a all

