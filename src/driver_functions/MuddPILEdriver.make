# MuddPILEdriver.make
# makes the master MuddPILEdriver program.
# make with: make -f MuddPILEdriver.make

CC = g++
CFLAGS= -c -I../../boost_mtl_minimal -Wall -O3 -std=c++11 -fPIC
OFLAGS = -I../../boost_mtl_minimal -Wall -O3 -std=c++11 -fPIC
LDFLAGS= -Wall -fPIC
SOURCES = MuddPILEdriver.cpp \
		../LSDRasterSpectral.cpp \
		../LSDIndexRaster.cpp \
		../LSDShapeTools.cpp \
		../LSDRaster.cpp \
		../LSDRasterModel.cpp \
		../LSDStatsTools.cpp \
		../LSDFlowInfo.cpp \
        ../LSDJunctionNetwork.cpp \
        ../LSDChannel.cpp \
        ../LSDIndexChannel.cpp \
        ../LSDMostLikelyPartitionsFinder.cpp \
		../LSDParticle.cpp \
        ../LSDRasterInfo.cpp \
        ../LSDSpatialCSVReader.cpp \
        ../LSDRasterMaker.cpp \
		../LSDLithoCube.cpp \
		../LSDParticleColumn.cpp \
        ../LSDParameterParser.cpp \
		../LSDCRNParameters.cpp
OBJ = $(SOURCES:.cpp=.o)
#LIBS = -lfftw3 -g -O0 -D_GLIBCXX_DEBUG
LIBS = -lfftw3 -Wwrite-strings
EXEC = MuddPILEdriver.out

all: $(SOURCES) $(SCRIPTS) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(OFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
