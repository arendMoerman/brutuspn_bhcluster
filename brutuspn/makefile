###################################################################
# Adjust path to mpreal.h file
MPRDIR  = /home/arend/Projects/Simulations/brutuspn_bhcluster/brutuspn/src/
###################################################################
LIB     = ./lib/
MPLIB   = ./mplib/
DIRSRC = ./src/

DIRTIM  = $(LIB)Timer/

DIRDH   = $(MPLIB)Data_Handler/

DIRRAN  = $(LIB)Random/
DIRTOO  = $(LIB)Tools/
DIRDIA  = $(MPLIB)Diagnostics/
DIRINIT = $(MPLIB)Initializer/

DIRPRO  = $(MPLIB)Products/
DIRDEL  = $(MPLIB)Delta/
###################################################################
COMP = g++
CFLAGS = -O3 -std=c++11
CC  = $(COMP) $(CFLAGS) -c
CO  = $(COMP) $(CFLAGS) -o
###################################################################
EXE = main.exe
all: $(EXE)
###################################################################
main.exe: main.o Brutus.o Bulirsch_Stoer.o Cluster.o Acceleration.o Star.o Delta.o Products.o Initializer.o Diagnostics.o Tools.o Random.o mtrand.o Data_Handler.o Timer.o  
	$(CO) main.exe main.o Brutus.o Bulirsch_Stoer.o Cluster.o Acceleration.o Star.o Delta.o Products.o Initializer.o Diagnostics.o Tools.o Random.o mtrand.o Data_Handler.o Timer.o  -lmpfr -lgmp
###################################################################
main.o: $(DIRSRC)main.cpp Brutus.o Bulirsch_Stoer.o Cluster.o Acceleration.o Star.o Delta.o Products.o Initializer.o Diagnostics.o Tools.o Random.o mtrand.o Data_Handler.o Timer.o  
	$(CC) $(DIRSRC)main.cpp -I$(MPRDIR) -I$(DIRTIM) -I$(DIRDH) -I$(DIRRAN) -I$(DIRDIA) -I$(DIRINIT) -I$(DIRPRO) -I$(DIRDEL) -I$(DIRTOO)
###################################################################
Brutus.o: $(DIRSRC)Brutus.h $(DIRSRC)Brutus.cpp Bulirsch_Stoer.o Cluster.o Acceleration.o Star.o Timer.o
	$(CC) $(DIRSRC)Brutus.cpp -I$(MPRDIR) -I$(DIRTIM) -I$(DIRSRC)

Bulirsch_Stoer.o: $(DIRSRC)Bulirsch_Stoer.h $(DIRSRC)Bulirsch_Stoer.cpp Cluster.o Acceleration.o Star.o Timer.o
	$(CC) $(DIRSRC)Bulirsch_Stoer.cpp -I$(MPRDIR) -I$(DIRTIM) -I$(DIRSRC)
Cluster.o: $(DIRSRC)Cluster.h $(DIRSRC)Cluster.cpp Acceleration.o Star.o 
	$(CC) $(DIRSRC)Cluster.cpp -I$(DIRINIT) -I$(MPRDIR) -I$(DIRSRC)
Acceleration.o: $(DIRSRC)Acceleration.h $(DIRSRC)Acceleration.cpp Star.o
	$(CC) $(DIRSRC)Acceleration.cpp -I$(MPRDIR) -I$(DIRSRC)
Star.o: $(DIRSRC)Star.h $(DIRSRC)Star.cpp 
	$(CC) $(DIRSRC)Star.cpp -I$(MPRDIR) -I$(DIRSRC)

Delta.o: $(DIRDEL)Delta.h $(DIRDEL)Delta.cpp
	$(CC) $(DIRDEL)Delta.cpp -I$(MPRDIR) -I$(DIRSRC)

Initializer.o: $(DIRINIT)Initializer.h $(DIRINIT)Initializer.cpp Data_Handler.o mtrand.o Random.o Tools.o Diagnostics.o Products.o
	$(CC) $(DIRINIT)Initializer.cpp -I$(DIRDH) -I$(DIRRAN) -I$(DIRDIA) -I$(DIRPRO) -I$(DIRTOO) -I$(MPRDIR)

Products.o: $(DIRPRO)Products.h $(DIRPRO)Products.cpp
	$(CC) $(DIRPRO)Products.cpp -I$(MPRDIR)
Diagnostics.o: $(DIRDIA)Diagnostics.h $(DIRDIA)Diagnostics.cpp Tools.o
	$(CC) $(DIRDIA)Diagnostics.cpp -I$(DIRTOO) -I$(MPRDIR)
Tools.o: $(DIRTOO)Tools.h $(DIRTOO)Tools.cpp
	$(CC) $(DIRTOO)Tools.cpp
Random.o: $(DIRRAN)Random.h $(DIRRAN)Random.cpp mtrand.o
	$(CC) $(DIRRAN)Random.cpp
mtrand.o: $(DIRRAN)mtrand.h $(DIRRAN)mtrand.cpp
	$(CC) $(DIRRAN)mtrand.cpp

Data_Handler.o: $(DIRDH)Data_Handler.h $(DIRDH)Data_Handler.cpp
	$(CC) $(DIRDH)Data_Handler.cpp -I$(MPRDIR)

Timer.o: $(DIRTIM)Timer.h $(DIRTIM)Timer.cpp
	$(CC) $(DIRTIM)Timer.cpp
###################################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe
	rm -f -r ccache
	rm -f .fuse*



