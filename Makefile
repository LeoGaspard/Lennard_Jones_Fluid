CC = g++
EXEC = Project
LIBS = 
FLAGS = -fopenmp -Wall -std=c++11

all: main.o
	$(CC) *.o -o $(EXEC) $(LIBS) $(FLAGS)
main.o : cdynamic.o
	$(CC) main.cpp -c $(FLAGS)
cdynamic.o : capplication.o
	$(CC) CDynamic.cpp -c $(FLAGS)
capplication.o: 
	$(CC) CApplication.cpp -c $(FLAGS)



clear:
	rm -f *o
mr_proper:
	rm -f *o $(EXEC)
