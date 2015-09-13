## sparskitJava

This project enables calling Fortran routines from Sparskit in Java.
It also includes basic Java classes for different matrix formats.
The /ipynb folder includes a IPython notebook that records and evaluates 

### Run it
1. Execute the makefile in /jni in order to build the interface to Fortran. The makefile might need to be modified according to the system's paths, compiler and operating system. It currently creates a .dylib (MacOS) and uses gfortran.
2. Running a Java program that makes use of the Sparskit class requires the JVM argument "-Djava.library.path=jni" so that the Fortran interface in /jni can be accessed.
3. In order to run the Java test, the data needs to be added to the /data folder. It currently contains a list of file names that were used for previous performance experiments. The files are taken from the [University of Florida Sparse Matrix Collection](http://www.cise.ufl.edu/research/sparse/matrices/index.html). See /data/files.txt for more details.
