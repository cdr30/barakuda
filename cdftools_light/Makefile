# Makefile for CDFTOOLS 

# ( make.macro is a link that points to the file macro.xxx where 
#   xxx is representative of your machine )
# !!  $Rev: 324 $
# !!  $Date: 2010-05-21 18:46:16 +0200 (Fri, 21 May 2010) $
# !!  $Id: Makefile 324 2010-05-21 16:46:16Z molines $
# !!--------------------------------------------------------------


include make.macro


CDFTOOLS=CDFTOOLS-2.1



EXEC = ./bin/cdfpsi.x \
       ./bin/cdfmaxmoc.x \
       ./bin/cdfmhst.x \
       ./bin/cdfsigtrp.x \
       ./bin/cdficediags.x \
       ./bin/cdfmoc.x \
       ./bin/cdfvT.x \
       ./bin/cdftransportiz.x

all: $(EXEC)

## Statistical programs
./bin/cdfmoy.x: src/cdfio.o   src/cdfmoy.f90
	$(F90) src/cdfmoy.f90 -o ./bin/cdfmoy.x  src/cdfio.o  $(FFLAGS)

./bin/cdfpsi.x: src/cdfio.o src/io_ezcdf.o src/cdfpsi.f90
	$(F90) src/cdfpsi.f90 -o ./bin/cdfpsi.x  src/cdfio.o src/io_ezcdf.o $(FFLAGS)

./bin/cdfeke.x: src/cdfio.o  src/cdfeke.f90
	$(F90) src/cdfeke.f90 -o ./bin/cdfeke.x  src/cdfio.o $(FFLAGS)

./bin/cdfrmsssh.x: src/cdfio.o  src/cdfrmsssh.f90
	$(F90) src/cdfrmsssh.f90 -o ./bin/cdfrmsssh.x  src/cdfio.o $(FFLAGS)

./bin/cdfvT.x: src/cdfio.o  src/io_ezcdf.o src/cdfvT.f90
	$(F90) src/cdfvT.f90 -o ./bin/cdfvT.x  src/cdfio.o src/io_ezcdf.o $(FFLAGS)

./bin/cdfcurl.x: src/cdfio.o  src/cdfcurl.f90
	$(F90) src/cdfcurl.f90 -o ./bin/cdfcurl.x  src/cdfio.o $(FFLAGS)

./bin/cdfmhst.x: src/cdfio.o  src/cdfmhst.f90
	$(F90) src/cdfmhst.f90 -o ./bin/cdfmhst.x  src/cdfio.o $(FFLAGS)

./bin/cdfmxl.x: src/cdfio.o src/eos.o  src/cdfmxl.f90
	$(F90) src/cdfmxl.f90 -o ./bin/cdfmxl.x  src/cdfio.o src/eos.o $(FFLAGS)

./bin/cdftransportiz.x: src/cdfio.o src/io_ezcdf.o src/cdftransportiz.f90
	$(F90) src/cdftransportiz.f90 -o ./bin/cdftransportiz.x  src/cdfio.o src/io_ezcdf.o $(FFLAGS)

./bin/cdficeflux.x: src/cdfio.o src/io_ezcdf.o src/cdficeflux.f90
	$(F90) src/cdficeflux.f90 -o ./bin/cdficeflux.x  src/cdfio.o $(FFLAGS)

./bin/cdfsigtrp.x: src/cdfio.o src/io_ezcdf.o src/eos.o  src/cdfsigtrp.f90
	$(F90)  src/cdfsigtrp.f90 -o ./bin/cdfsigtrp.x  src/cdfio.o src/io_ezcdf.o src/eos.o $(FFLAGS)

./bin/cdfmoc.x: src/cdfio.o src/io_ezcdf.o src/cdfmoc.f90
	$(F90) src/cdfmoc.f90 -o ./bin/cdfmoc.x  src/cdfio.o src/io_ezcdf.o $(FFLAGS)

./bin/cdfmocsig.x: src/cdfio.o src/eos.o  src/cdfmocsig.f90
	$(F90) src/cdfmocsig.f90 -o ./bin/cdfmocsig.x  src/cdfio.o src/eos.o $(FFLAGS)

./bin/cdfmocatl.x: src/cdfio.o  src/cdfmocatl.f90
	$(F90) src/cdfmocatl.f90 -o ./bin/cdfmocatl.x  src/cdfio.o $(FFLAGS)

./bin/cdfmean.x: src/cdfio.o  src/cdfmean.f90
	$(F90) src/cdfmean.f90 -o ./bin/cdfmean.x  src/cdfio.o $(FFLAGS)

./bin/cdfmeanvar.x: src/cdfio.o  src/cdfmeanvar.f90
	$(F90) src/cdfmeanvar.f90 -o ./bin/cdfmeanvar.x  src/cdfio.o $(FFLAGS)

./bin/cdfheatc.x: src/cdfio.o  src/cdfheatc.f90
	$(F90) src/cdfheatc.f90 -o ./bin/cdfheatc.x  src/cdfio.o $(FFLAGS)

./bin/cdfmxlheatc.x: src/cdfio.o  src/cdfmxlheatc.f90
	$(F90) src/cdfmxlheatc.f90 -o ./bin/cdfmxlheatc.x  src/cdfio.o $(FFLAGS)

./bin/cdficediags.x: src/cdfio.o  src/cdficediags.f90
	$(F90) src/cdficediags.f90 -o ./bin/cdficediags.x  src/cdfio.o $(FFLAGS)

./bin/cdfzonalmean.x: src/cdfio.o  src/cdfzonalmean.f90
	$(F90) src/cdfzonalmean.f90 -o ./bin/cdfzonalmean.x  src/cdfio.o $(FFLAGS) 

./bin/cdfzonalsum.x: src/cdfio.o  src/cdfzonalsum.f90
	$(F90) src/cdfzonalsum.f90 -o ./bin/cdfzonalsum.x  src/cdfio.o $(FFLAGS) 

./bin/cdfzonalout.x: src/cdfio.o  src/cdfzonalout.f90
	$(F90) src/cdfzonalout.f90 -o ./bin/cdfzonalout.x  src/cdfio.o $(FFLAGS) 

./bin/cdfhflx.x: src/cdfio.o  src/cdfhflx.f90
	$(F90) src/cdfhflx.f90 -o ./bin/cdfhflx.x  src/cdfio.o $(FFLAGS)

./bin/cdfbuoyflx: src/cdfio.o  src/eos.o cdfbuoyflx.f90
	$(F90) src/cdfbuoyflx.f90 -o ./bin/cdfbuoyflx.x  src/cdfio.o src/eos.o $(FFLAGS)

./bin/cdfmaxmoc.x: src/cdfio.o src/io_ezcdf.o src/cdfmaxmoc.f90
	$(F90) src/cdfmaxmoc.f90 -o ./bin/cdfmaxmoc.x  src/cdfio.o src/io_ezcdf.o $(FFLAGS)

./bin/cdfw.x: src/cdfio.o  src/cdfw.f90
	$(F90) src/cdfw.f90 -o ./bin/cdfw.x  src/cdfio.o $(FFLAGS)


## Modules

src/cdfio.o: src/cdfio.f90
	@mkdir -p mod
	@mkdir -p bin
	$(F90) -c src/cdfio.f90 $(FFLAGS) -o src/cdfio.o


src/io_ezcdf.o: src/io_ezcdf.f90
	@mkdir -p mod
	@mkdir -p bin
	$(F90) -c src/io_ezcdf.f90 $(FFLAGS) -o src/io_ezcdf.o



#src/modcdfnames.o: src/modcdfnames.f90
#	@mkdir -p mod
#	$(F90) -c src/modcdfnames.f90 $(FFLAGS) -o src/modcdfnames.o

src/eos.o: src/eos.f90
	@mkdir -p mod
	@mkdir -p bin
	$(F90) -c src/eos.f90   $(FFLAGS) -o src/eos.o

cdftools.o: src/cdfio.o cdftools.f90
	@mkdir -p mod
	@mkdir -p bin
	$(F90) -c cdftools.f90 $(FFLAGS)

modpoly.o: modpoly.f90
	@mkdir -p mod
	$(F90) -c modpoly.f90 $(FFLAGS)

## Utilities
tar:
	( cd ../ ; tar cf cdftools-2.1.tar $(CDFTOOLS)/*90 $(CDFTOOLS)/Make* \
          $(CDFTOOLS)/section.dat $(CDFTOOLS)/JOBS $(CDFTOOLS)/DOC \
          $(CDFTOOLS)/macro.* )

clean:
	\rm -rf mod src/*.o  *~ bin

cleanexe: clean
	\rm -f $(EXEC)

#install:
#	\mv -f $(EXEC) ../bin/
#