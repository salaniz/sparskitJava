# Define a variable for classpath
CLASS_PATH = ../target/classes/

#Define package path
PCK_PATH = edu/tuberlin/sparskit

# Define package name
PCK_NAME = edu.tuberlin.sparskit.

# Define a virtual path for .class in the bin directory
vpath %.class $(CLASS_PATH)$(PCK_PATH)

FC := gfortran

# On windows reanme to dll
# OUTPUTNAME=sparskit.dll
OUTPUTNAME=libsparskit.so

#PLATFORM_INCLUDE=${JAVA_HOME}/include/win32
PLATFORM_INCLUDE=${JAVA_HOME}/include/linux

all : ${OUTPUTNAME}

# $@ matches the target, $^ matches all dependancies
${OUTPUTNAME} : SparskitInterface.o SparskitMatVecFortran.o SparskitFormatFortran.o SparskitUnaryFortran.o SparskitBlassmFortran.o SparskitInfofunFortran.o
	gcc -m64 -shared $^ -o $@ -lgfortran

# $@ matches the target, $< matches the first dependancy
SparskitMatVecFortran.o : matvec.f
	$(FC) -m64 -fPIC -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitFormatFortran.o : formats.f
	$(FC) -m64 -fPIC -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitUnaryFortran.o : unary.f
	$(FC) -m64 -fPIC -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitBlassmFortran.o : blassm.f
	$(FC) -m64 -fPIC -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitInfofunFortran.o : infofun.f
	$(FC) -m64 -fPIC -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitInterface.o : SparskitInterface.c edu_tuberlin_sparskit_Sparskit.h
	gcc -m64 -fPIC -std=c99 -I"${JAVA_HOME}/include"  -I"${PLATFORM_INCLUDE}" -c $< -o $@

# $* matches the target filename without the extension
edu_tuberlin_sparskit_Sparskit.h : Sparskit.class
	${JAVA_HOME}/bin/javah -classpath $(CLASS_PATH) $(PCK_NAME)Sparskit

clean :
	rm -f *.o *.h *.so
