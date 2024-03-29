## sparskitJava

This project enables calling Fortran routines from Sparskit in Java.
It also includes basic Java classes for different matrix formats.
The /ipynb folder includes a IPython notebook that records and evaluates performance on a set of sample matrices.

### Requirements
* Java 1.8
* gfortran
* Maven

## Get the sample Files
1. execute `download.sh` in the data folder

> Note for running on Windows

Get TDM64 http://tdm-gcc.tdragon.net/download
 - install with fortran for your architecture
 - compile and rename resulting file to
 - restart the IDE

### Run it
1. Install Maven dependencies in pom.xml and make sure the Java files are compiled.
   The class files should be located in ./target/classes.
2. Execute the makefile in /jni in order to build the interface to Fortran.
   The makefile might need to be modified according to the system's paths, compiler and operating system.
   The main makefile runs on Linux and uses gfortran. There is another makefile for MacOS users (makefile_macos).
3. Running a Java program that makes use of the Sparskit class requires the JVM argument "-Djava.library.path=jni"
   so that the Fortran interface in /jni can be accessed.
4. In order to run the Java test, the data needs to be added to the /data folder.
   It currently contains a list of file names that were used for previous performance experiments.
   The files are taken from the [University of Florida Sparse Matrix Collection](http://www.cise.ufl.edu/research/sparse/matrices/index.html).
   See /data/files.txt for more details.

### Open work
The file futurework.txt contains a list of known issues and open work/questions about the project.
