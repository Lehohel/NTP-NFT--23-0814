# IDIR=./
CXX = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.5\bin\nvcc.exe"

CXXFLAGS += $(shell pkg-config --cflags --libs opencv4)
LDFLAGS += $(shell pkg-config --libs --static opencv)
FILENAME = cuABM

all: clean build run

build:
	$(CXX) .\src\$(FILENAME).cu -o .\bin\$(FILENAME).exe -I..\Common -I.\

run:
	.\bin\$(FILENAME).exe $(ARGS)

clean:
	del ".\bin\$(FILENAME).exe"
	del ".\bin\$(FILENAME).exp"
	del ".\bin\$(FILENAME).lib"