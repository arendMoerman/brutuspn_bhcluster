###################################################################
ROOT    = /data1/boekholt/2sirius/

LIB     = $(ROOT)lib/

DIRDH   = $(LIB)Data_Handler/
DIRRAN  = $(LIB)Random/
DIRDIA  = $(LIB)Diagnostics/
DIRPRO  = $(LIB)Products/
DIRINIT = $(LIB)Initializer/
###################################################################
CC  = g++ -O3 -c
CO  = g++ -O3 -o
###################################################################
EXE = test.exe
all: $(EXE)
###################################################################
test.exe: test.o Data_Handler.o mtrand.o Random.o Diagnostics.o Products.o Initializer.o
	$(CO) test.exe test.o Data_Handler.o mtrand.o Random.o Diagnostics.o Products.o Initializer.o
###################################################################
test.o: test.cpp Data_Handler.o mtrand.o Random.o Diagnostics.o Products.o Initializer.o
	$(CC) test.cpp -I$(DIRDH) -I$(DIRRAN) -I$(DIRDIA) -I$(DIRINIT) -I$(DIRPRO)

Initializer.o: $(DIRINIT)Initializer.h $(DIRINIT)Initializer.cpp Data_Handler.o mtrand.o Random.o Diagnostics.o Products.o
	$(CC) $(DIRINIT)Initializer.cpp -I$(DIRDH) -I$(DIRRAN) -I$(DIRDIA) -I$(DIRPRO)
Products.o: $(DIRPRO)Products.h $(DIRPRO)Products.cpp
	$(CC) $(DIRPRO)Products.cpp
Diagnostics.o: $(DIRDIA)Diagnostics.h $(DIRDIA)Diagnostics.cpp
	$(CC) $(DIRDIA)Diagnostics.cpp
Random.o: $(DIRRAN)Random.h $(DIRRAN)Random.cpp mtrand.o
	$(CC) $(DIRRAN)Random.cpp
mtrand.o: $(DIRRAN)mtrand.h $(DIRRAN)mtrand.cpp
	$(CC) $(DIRRAN)mtrand.cpp
Data_Handler.o: $(DIRDH)Data_Handler.h $(DIRDH)Data_Handler.cpp
	$(CC) $(DIRDH)Data_Handler.cpp
###################################################################
clean:
	rm -f *~
	rm -f *.o
	rm -f *.exe
	rm -f -r ccache
	rm -f .fuse*



