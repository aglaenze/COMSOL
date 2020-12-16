INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Install/lib

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield

all : plotField spectrumFe55 avalanche signal

plotField: plotField.C 
	$(CXX) $(CFLAGS) plotField.C
	$(CXX) -o plotField plotField.o $(LDFLAGS)
	rm plotField.o

spectrumFe55: spectrumFe55.C 
	$(CXX) $(CFLAGS) spectrumFe55.C
	$(CXX) -o spectrumFe55 spectrumFe55.o $(LDFLAGS)
	rm spectrumFe55.o

avalanche: avalanche.C 
	$(CXX) $(CFLAGS) avalanche.C
	$(CXX) -o avalanche avalanche.o $(LDFLAGS)
	rm avalanche.o

signal: signal.C 
	$(CXX) $(CFLAGS) signal.C
	$(CXX) -o signal signal.o $(LDFLAGS)
	rm signal.o


