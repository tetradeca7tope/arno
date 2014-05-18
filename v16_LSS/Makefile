JAVA_HOME=/usr/lib/jvm/java-6-openjdk-amd64/
INC=-I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux

all: liblss_wrapper.so LssWrapper.class

LssWrapper.class: LssWrapper.java
	javac LssWrapper.java

LssWrapper.h: LssWrapper.class LssWrapper.java
	javah LssWrapper

liblss_wrapper.so: LssWrapper.o compute_lrg_likelihood.o
	gcc -shared LssWrapper.o compute_lrg_likelihood.o -lg2c -o liblss_wrapper.so

LssWrapper.o: LssWrapper.c LssWrapper.h
	gcc $(INC) -c LssWrapper.c


compute_lrg_likelihood.o: compute_lrg_likelihood.f
	gfortran -c compute_lrg_likelihood.f

a.out: compute_lrg_likelihood.f
	gfortran compute_lrg_likelihood.f

clean:  
	rm -f *.o *.class a.out *.so *~

