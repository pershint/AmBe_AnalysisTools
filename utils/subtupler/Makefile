ROOT_CONFIG := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
RATROOT := $(shell echo $(RATROOT))
RAT_LIBRARY_PATH := $(RATROOT)/lib
RAT_INCLUDE_PATH := $(RATROOT)/include
RAT_LIBRARY := RATEvent_Linux
PROGRAM := subtupleMaker_AmBe.cpp
OUT := subtupleMaker_AmBe

standard:
	g++ -fpermissive -std=c++11 $(ROOT_CONFIG) -I$(RAT_INCLUDE_PATH) -I$(RAT_INCLUDE_PATH)/libpq -L$(RAT_LIBRARY_PATH) $(PROGRAM) -o $(OUT) -l$(RAT_LIBRARY) $(ROOT_LIBS)
