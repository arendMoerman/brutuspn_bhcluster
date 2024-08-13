###################################################################
HOME    = /data1/boekholt/

DIRMPI  = $(HOME)Packages/mpi_local/bin/

MPRDIR  = $(HOME)Packages/mpreal_local/include/
MPRLIB  = $(HOME)Packages/mpreal_local/lib/

ROOT    = $(HOME)2sirius/

LIB     = $(ROOT)lib/
MPLIB   = $(ROOT)mplib/

DIRCOM  = $(MPLIB)mpCommunicator/
###################################################################
COMP = $(DIRMPI)mpicxx
CFLAGS = -O1
CC  = $(COMP) $(CFLAGS) -c
CO  = $(COMP) $(CFLAGS) -o
###################################################################
EXE = test.exe
all: $(EXE)
###################################################################
test.exe: test.o mpCommunicator.o 
	$(CO) test.exe test.o mpCommunicator.o $(MPRLIB)libmpfr.a $(MPRLIB)libgmp.a
###################################################################
test.o: test.cpp mpCommunicator.o 
	$(CC) test.cpp -I$(MPRDIR) -I$(DIRCOM)
mpCommunicator.o: $(DIRCOM)mpCommunicator.h $(DIRCOM)mpCommunicator.cpp 
	$(CC) $(DIRCOM)mpCommunicator.cpp -I$(MPRDIR) -I$(DIRMPI)
###################################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe
	rm -f -r ccache
	rm -f .fuse*



