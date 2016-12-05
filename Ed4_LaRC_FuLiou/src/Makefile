PROG = fu_interface

MOD_FLAG = -I../lib/mod 
LIBS = ../lib/libEd3Fu_201212.a
FUVER=Ed4_LaRC_FuLiou/src

PYVERSION:=35m
#set python to version 3.5
PLATFORM:=$(shell uname -s)
 
ifeq ($(PLATFORM),Darwin)
	PYLIB:=fupy.cpython-$(PYVERSION)-$(shell echo $(PLATFORM) | tr A-Z a-z).so
endif
ifeq ($(PLATFORM),Linux)
	PYLIB:=fupy.cpython-$(PYVERSION)-x86_64-$(shell echo $(PLATFORM) | tr A-Z a-z)-gnu.so
endif

 

#CC = 
#CFLAGS = 

#!!!SGI
F90 =gfortran -fdefault-real-8 -fdefault-double-8 -fbacktrace -fno-align-commons
F77 =gfortran -fdefault-real-8 -fdefault-double-8 -fbacktrace -fno-align-commons
F90LOAD =
F90COMP = -c  
FCOMP	= -c  

#!!!MAC IBM XLF
#F90 = f90  
#F77 = -qxlf90=nosignedzero,autodealloc -O2 -qmaxmem=65536 -c -qextname -qsuffix=f=f90 -q64 -qwarn64 -qinit=f90ptr
#F90COMP = -c -qxlf90=nosignedzero,autodealloc -O2 -qmaxmem=65536 -c -qextname -qsuffix=f=f90 -q64 -qwarn64 -qinit=f90ptr
#FCOMP   = -c -qxlf90=nosignedzero,autodealloc -O2 -qmaxmem=65536 -c -qextname -qsuffix=f=f90 -q64 -qwarn64 -qinit=f90ptr
#F90LOAD = 

# F90 = xlf90  -pg
# F90COMP = -c -qsuffix=f=f90  
# FCOMP = -c -qextname
# F90LOAD= -static

all: $(PROG)
	f2py -c --f90flags="-Wno-tabs" -m fupy fupy.pyf fu_interface.f90 $(LIBS) $(MOD_FLAG) 
	install_name_tool $(PYLIB) -change $(LIBS) $(realpath $(LIBS))
	otool -L $(PYLIB)
	cd ../; ln -sf $(FUVER)/$(PYLIB)



$(PROG):   $(PROG).o $(LIBS) $(CLIBS) 
	 $(F90) $(F90LOAD) -o $@ $@.o $(LIBS) $(CLIBS) 


clean:
	-\rm -f $(PROG) 
	-\rm -f *.o
	-\rm -f *.mod 
	

.SUFFIXES: 

.SUFFIXES: .f90 .mod .f .c .o

.f90.o:
	$(F90) $(F90COMP) $(MOD_FLAG) $(INC_FLAG) $< 

.f90.mod:
	$(F90) $(F90COMP) $(MOD_FLAG) $(INC_FLAG) $< 

.f.o:
	$(F90) $(FCOMP)  $(MOD_FLAG) $(INC_FLAG) $< 

.c.o:
	$(CC) $(CFLAGS)  $(MOD_FLAG) $(INC_FLAG) $< 


simple : simple.o

