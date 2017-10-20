# Makefile for YUM

#FC    = ifort
FC    = gfortran
FFLAGS    = -O3
SRC    = yum-dp-orc2.f90
EXE    = a.out

###INCLUDE = -I /opt/netcdf/include
INCLUDE =
###LIB    = /opt/netcdf/lib/libnetcdf.a
LIB    =
LOG = ../log
TIMECMD    = time
COMPILELOG    = $(LOG)/compile.`date '+%Y%m%d.%H%M%S'`.log
RUNLOG    = $(LOG)/run.`date '+%Y%m%d.%H%M%S'`.log
TARFILE = yumsrc.`date '+%Y%m%d.%H%M%S'`.tgz

$(EXE)    : $(SRC) $(HEADER)
	($(FC) $(FFLAGS) $(SRC) -o $(EXE) $(INCLUDE) $(LIB)) 2>&1 | tee $(COMPILELOG)

run    :
	$(TIMECMD) ./$(EXE) 2>&1 | tee $(RUNLOG)

clean    :
	rm -f $(EXE) core

tar    :
	tar zcvf $(TARFILE) $(SRC) $(HEADER) Makefile

backup    :
	cp $(SRC) .$(SRC).`date '+%Y%m%d%H%M%S'`
	cp Makefile .Makefile.`date '+%Y%m%d%H%M%S'`

xxx    :
	cp /dev/null log.`date '+%Y%m%d%H%M%S'`
