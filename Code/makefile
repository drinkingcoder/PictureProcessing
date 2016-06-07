main:main.o bmpFile.o Matrix.o
	g++ main.o bmpFile.o Matrix.o -o main
	rm *.o

main.o:main.cpp bmpFile.h
	g++ main.cpp -c

bmpFile.o:bmpFile.cpp bmpFile.h
	g++ bmpFile.cpp -c

Matrix.o:Matrix/Matrix.h Matrix/Matrix.cpp
	g++ Matrix/Matrix.cpp -c 

clean:
	rm *.o
	rm main
	rm Matrix/Matrix.o
