# Define a variable for classpath
CLASS_PATH = ../target/classes/

#Define package path
PCK_PATH = edu/tuberlin/sparskit

# Define package name
PCK_NAME = edu.tuberlin.sparskit.

# Define a virtual path for .class in the bin directory
vpath %.class $(CLASS_PATH)$(PCK_PATH)

FC := /usr/local/bin/gfortran

all : libsparskit.dylib

# $@ matches the target, $^ matches all dependancies
libsparskit.dylib : SparskitInterface.o SparskitMatVecFortran.o SparskitFormatFortran.o SparskitUnaryFortran.o SparskitBlassmFortran.o SparskitInfofunFortran.o
	gcc -m64 -L"/usr/local/gfortran/lib" -Wl -shared $^ -o $@ -lgfortran

# $@ matches the target, $< matches the first dependancy
SparskitMatVecFortran.o : matvec.f
	$(FC) -m64 -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitFormatFortran.o : formats.f
	$(FC) -m64 -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitUnaryFortran.o : unary.f
	$(FC) -m64 -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitBlassmFortran.o : blassm.f
	$(FC) -m64 -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitInfofunFortran.o : infofun.f
	$(FC) -m64 -c $< -o $@
	
# $@ matches the target, $< matches the first dependancy
SparskitInterface.o : SparskitInterface.c edu_tuberlin_sparskit_Sparskit.h
	gcc -m64 -I"/Library/Java/JavaVirtualMachines/jdk1.8.0_25.jdk/Contents/Home/include" -c $< -o $@

# $* matches the target filename without the extension
edu_tuberlin_sparskit_Sparskit.h : Sparskit.class
	javah -classpath $(CLASS_PATH) $(PCK_NAME)Sparskit

clean :
	rm *.o *.h *.dylib