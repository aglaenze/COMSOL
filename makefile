INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield


LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm
# LDFLAGS += -g

all : plotField spectrumFe55 drawGeometry gain avalanche signal ibf

plotField: plotField.C 
	$(CXX) $(CFLAGS) plotField.C
	$(CXX) `root-config --cflags` -o plotField plotField.o $(LDFLAGS)
	rm plotField.o

spectrumFe55: spectrumFe55.C 
	$(CXX) $(CFLAGS) spectrumFe55.C
	$(CXX) `root-config --cflags` -o spectrumFe55 spectrumFe55.o $(LDFLAGS)
	rm spectrumFe55.o

drawGeometry: drawGeometry.C 
	$(CXX) $(CFLAGS) drawGeometry.C
	$(CXX) `root-config --cflags` -o drawGeometry drawGeometry.o $(LDFLAGS)
	rm drawGeometry.o

gain: gain.C 
	$(CXX) $(CFLAGS) gain.C
	$(CXX) `root-config --cflags` -o gain gain.o $(LDFLAGS)
	rm gain.o

avalanche: avalanche.C 
	$(CXX) $(CFLAGS) avalanche.C
	$(CXX) `root-config --cflags` -o avalanche avalanche.o $(LDFLAGS)
	rm avalanche.o

signal: signal.C 
	$(CXX) $(CFLAGS) signal.C
	$(CXX) `root-config --cflags` -o signal signal.o $(LDFLAGS)
	rm signal.o

ibf: ibf.C 
	$(CXX) $(CFLAGS) ibf.C
	$(CXX) `root-config --cflags` -o ibf ibf.o $(LDFLAGS)
	rm ibf.o

