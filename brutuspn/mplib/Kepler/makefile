###################################################################
ROOT    = /data1/boekholt/

MPRDIR  = $(ROOT)Packages/mpreal_local/include/
MPRLIB  = $(ROOT)Packages/mpreal_local/lib/

HOME    = $(ROOT)2sirius/

LIB     = $(HOME)lib/
MPLIB   = $(HOME)mplib/

DIRKEP  = $(MPLIB)Kepler/
###################################################################
COMP = g++
CFLAGS = -O3
CC  = $(COMP) $(CFLAGS) -c
CO  = $(COMP) $(CFLAGS) -o
###################################################################
EXE = test.exe test_nbody.exe
all: $(EXE)
###################################################################
test.exe: test.o Kepler.o   
	$(CO) test.exe test.o Kepler.o $(MPRLIB)libmpfr.a $(MPRLIB)libgmp.a
test_nbody.exe: test_nbody.o Kepler.o   
	$(CO) test_nbody.exe test_nbody.o Kepler.o $(MPRLIB)libmpfr.a $(MPRLIB)libgmp.a
###################################################################
test.o: test.cpp Kepler.o  
	$(CC) test.cpp -I$(MPRDIR) -I$(DIRKEP)
test_nbody.o: test_nbody.cpp Kepler.o  
	$(CC) test_nbody.cpp -I$(MPRDIR) -I$(DIRKEP)
###################################################################
Kepler.o: $(DIRKEP)Kepler.h $(DIRKEP)Kepler.cpp
	$(CC) $(DIRKEP)Kepler.cpp -I$(MPRDIR)
###################################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe
	rm -f -r ccache
	rm -f .fuse*



