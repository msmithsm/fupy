PROG = 

#SRCS =  extras.f90 fuinput.f90 fuoutput.f90 fuprint.f90 \
#	icedirsfc.f90 ma_tip.f90 qsortd.f90 vla.f90 rad_multi_200511.f90 seiji_k2d.f90 \
#	seiji_solver_200511.f90  taucorr.f90 uvcor_all.f90 \
#	 gflq.f90  calipso_output.f90 zjin.f90 ar_asy.f90

OBJS =  extras.o fuinput.o fuoutput.o fuprint.o entropy_lw.o icedirsfc.o ma_tip.o \
	qsortd.o vla.o rad_multi_200511.o seiji_k2d.o seiji_solver_200511.o  \
	taucorr.o uvcor_all.o \
	gflq.o calipso_output.o zjin.o ar_asy.o

OBJS_F = aqua_wnflt_0404.o   misc_200511.o \
	cloud_optics.o 	seiji_twostreamsolv_sw_v20.o aerosols_200511.o

##X86
#OBJS_C = qsort_c.o

FULIB = libEd3Fu_201212.a 
LIBS = 
MOD_FLAG = 
INC_FLAG = 

#CC = 
#CFLAGS = 

#!!!SGI
F90 = gfortran -O3 -fdefault-real-8 -fdefault-double-8 -fno-align-commons  
F77 = gfortran -O3 -fdefault-real-8 -fdefault-double-8 -fno-align-commons 
F90LOAD = 
F90COMP = -c  
FCOMP   = -c   

#!!!MAC IBM XLF
#F90 = f90  -qextname=qsort
#F77 = xlf  
## F90COMP = -pg  -c
## FCOMP   = -pg  -c
#F90LOAD =  
#FCOMP=-O2 -c -qmaxmem=8192 -qextname 

#  F90 = xlf90  -pg
# F90COMP = -c -qsuffix=f=f90  -qextname
# FCOMP = -c -qextname

# THESE ARE NONSTANDARD for PPC...
#F90COMP = -qxlf90=nosignedzero,autodealloc -O5 -lmass -qmaxmem=65536 -c -qextname -qsuffix=f=f90 -q64 -qwarn64 -qinit=f90ptr 
##FCOMP = -qxlf90=nosignedzero,autodealloc -O5 -lmass -qmaxmem=32768 -c -qextname -qfixed=80 -q64 -qwarn64
#F90LOAD = -q64 
#

all: $(FULIB) 

$(FULIB): $(OBJS) $(OBJS_F)  $(OBJS_C)
	 ar -rcv $(FULIB) $(OBJS) $(OBJS_F) $(OBJS_C) 
	-\cp *.mod ../mod 
	-\cp $(FULIB) ../ 
tidy:
	-\rm -f $(PROG) $(OBJS) $(OBJS_F) $(OBJS_C)
	-\rm -f *.mod 
	-\rm -f *.a 
clean:
	-\rm -f $(PROG) $(OBJS) $(OBJS_F) $(OBJS_C)
	-\rm -f *.mod 
	-\rm -f *.mod 
	-\rm -f *.a 
	-\rm -f ../mod/*.mod 
	-\rm -f ../$(FULIB)

.SUFFIXES: 

.SUFFIXES: .f90 .mod .f .c .o

.f90.o:
	$(F90) $(F90COMP)  $(MOD_FLAG) $(INC_FLAG) $< 

.f90.mod:
	$(F90) $(F90COMP)  $(MOD_FLAG) $(INC_FLAG) $< 

.f.o:
	$(F77) $(FCOMP)  $(MOD_FLAG) $(INC_FLAG) $< 

.c.o:
	$(CC) $(CFLAGS)  $(MOD_FLAG) $(INC_FLAG) $< 



fuinput.o: taucorr.o extras.o
fuoutput.o: fuinput.o
fuprint.o: fuinput.o fuoutput.o
icedirsfc.o: fuinput.o fuoutput.o
ma_tip.o: fuinput.o
rad_multi_200511.o: fuinput.o fuoutput.o icedirsfc.o uvcor_all.o
seiji_k2d.o: fuinput.o 
seiji_solver_200511.o: fuinput.o fuoutput.o
uvcor_all.o: fuinput.o fuoutput.o 
wpar.o: rad_multi_200511.o
cloud_optics.o : fuinput.o

vla.o : qsortd.o
