# IDIR=./
CXX = cl

CXXFLAGS += $(shell pkg-config --cflags --libs opencv4)
LDFLAGS += $(shell pkg-config --libs --static opencv)
FILENAME = stABM

all: clean build run

build:
	$(CXX) .\src\$(FILENAME).cpp -Fe.\bin\$(FILENAME).exe -I..\Common -I.\ -EHsc -O2

run:
	.\bin\$(FILENAME).exe $(ARGS)

clean:
	del ".\bin\$(FILENAME).exe"