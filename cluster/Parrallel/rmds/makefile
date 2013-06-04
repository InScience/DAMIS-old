vpath %.cpp ../libs
CC=mpiCC
CPPFLAGS=-I /home/Regimantas/cpp -I /home/Regimantas/cpp/libs -I /home/Regimantas/cpp/libs/alglib

PROJ_BIN=rmds.out
PROJ_SRC=rmdsmain.cpp
OBJECTDIR=build

all:build ${PROJ_BIN}

build:

# target : prerequisite (butinos salygos)
#
${PROJ_BIN}: ${OBJECTDIR}/mds.o ${OBJECTDIR}/rmds.o ${OBJECTDIR}/lbfgs.o ${OBJECTDIR}/alglib_pca.o ${OBJECTDIR}/tools.o \
${OBJECTDIR}/ap.o ${OBJECTDIR}/correlation.o ${OBJECTDIR}/pca.o ${OBJECTDIR}/blas.o ${OBJECTDIR}/descriptivestatistics.o ${OBJECTDIR}/svd.o \
${OBJECTDIR}/sblas.o ${OBJECTDIR}/ialglib.o ${OBJECTDIR}/hblas.o ${OBJECTDIR}/reflections.o ${OBJECTDIR}/creflections.o \
${OBJECTDIR}/ablasf.o ${OBJECTDIR}/ablas.o ${OBJECTDIR}/ortfac.o ${OBJECTDIR}/rotations.o ${OBJECTDIR}/bdsvd.o

# -MMD -MP -MF $@.d áraðo priklausomybes i build.d failà
# comanda turi buti atskiriama tabais
# rankytemis is to failo isimame headeris nuo kuriu priklauso kad sukurti reikiamus objektus
	${CC} -O2 -MMD -MP -MF $@.d ${CPPFLAGS} -o ${PROJ_BIN} $^ ${PROJ_SRC}
#	${CC} ${CPPFLAGS} -o ${PROJ_BIN} $^ ${PROJ_SRC}


${OBJECTDIR}/mds.o: ../mds.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/rmds.o: ../rmds.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/dma.o: ../dma.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/lbfgs.o: ../lbfgs.cpp ../lbfgs.h
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/alglib_pca.o: ../alglib_pca.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/tools.o: ../tools.cpp tools.d
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/%.o: ../libs/%.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

${OBJECTDIR}/%.o: ../libs/alglib/%.cpp
	${CC} ${CPPFLAGS} -c $< -o $@

%.d: ../%.cpp
	$(CC) $(CPPFLAGS) -M $< | sed 's/\($*\)\.o[ :]*/\1.o $@: /g' > $@

libs/%.d: ../libs/%.cpp
	$(CC) $(CPPFLAGS) -M $< | sed 's/\($*\)\.o[ :]*/\1.o $@: /g' > $@
	


# arba
#	$(CC) -M $(CPPFLAGS) $< > $@.$$$$
#	sed 's/\($*\)\.o[ :]*,\1.o $@ : /g' < $@.$$$$ > $@
#	rm -f $@.$$$$


${OBJECTDIR}/ap2.o:
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(CC) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ap.o ../libs/alglib/ap.cpp


${OBJECTDIR}/tools2.o: ../tools.h  ../libs/alglib/ap.h
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(CC) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/tools.o ../tools.cpp ../libs/alglib/ap.cpp ../libs/alglib/correlation.cpp

clean:
	rm -f *.core ${PROJ_BIN}
	rm -f ${OBJECTDIR}/*.o 

rebuild: clean build ${PROJ_BIN}
	
