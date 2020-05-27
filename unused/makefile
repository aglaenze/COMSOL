INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Install/lib

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield

all : drawGeometry

drawGeometry: drawGeometry.C 
	$(CXX) $(CFLAGS) drawGeometry.C
	$(CXX) -o drawGeometry drawGeometry.o $(LDFLAGS)
	rm drawGeometry.o




