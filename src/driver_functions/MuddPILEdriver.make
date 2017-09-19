# MuddPILEdriver.make 
# makes the master MuddPILEdriver program. 
# make with: make -f MuddPILEdriver.make

CC = g++
CFLAGS= -c -I../../boost_mtl_minimal -Wall -O3
OFLAGS = -I../../boost_mtl_minimal -Wall -O3
LDFLAGS= -Wall
SOURCES = MuddPILEdriver.cpp \
		../LSDRasterSpectral.cpp \
		../LSDIndexRaster.cpp \
		../LSDShapeTools.cpp \
		../LSDRaster.cpp \
		../LSDRasterModel.cpp \
		../LSDStatsTools.cpp \
		../LSDFlowInfo.cpp \
		../LSDParticle.cpp \
        ../LSDRasterMaker.cpp \
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
