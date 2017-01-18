#DATE = `date "+%Y%m%d.%H%M"`
DATE = `date "+%Y%m%d"`
NAME = "Ed4_LaRC_FuLiou"
all: fuinterface

libsrc:
	cd ./lib/src ;make

simple: libsrc
	cd ./src/simple ; make
	cd ./src/simple ; simple

fuinterface : libsrc
	cd ./src/ ; make

clean : 
	cd ./lib/src ;make clean
	cd ./src/simple ; make clean
	cd ./src ; make clean
	rm -f *.so *.a *.o 

#------------------------------------
tar:
	echo ${DATE}
	tar cvf ../${NAME}${DATE}.tar . 
	gzip    ../${NAME}${DATE}.tar
	ls -l   ../${NAME}${DATE}.tar.gz
