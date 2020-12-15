CC = g++
EXEC = Project
LIBS = 
FLAGS = -fopenmp -Wall -std=c++2a  -O3

all: main.o
	$(CC) *.o -o $(EXEC) $(LIBS) $(FLAGS)
main.o : cdynamic.o
	$(CC) src/main.cpp -c $(FLAGS) $(LIBS)
cdynamic.o : capplication.o
	$(CC) src/CDynamic.cpp -c $(FLAGS) $(LIBS)
capplication.o: cbox.o
	$(CC) src/CApplication.cpp -c $(FLAGS) $(LIBS)
cbox.o: catom.o c3mat.o
	$(CC) src/CBox.cpp -c $(FLAGS) $(LIBS)
catom.o: cpos.o cspeed.o cforce.o
	$(CC) src/CAtom.cpp -c $(FLAGS) $(LIBS)
cpos.o: c3vec.o
	$(CC) src/CPos.cpp -c $(FLAGS) $(LIBS)
cspeed.o: c3vec.o
	$(CC) src/CSpeed.cpp -c $(FLAGS) $(LIBS)
cforce.o: c3vec.o
	$(CC) src/CForce.cpp -c $(FLAGS) $(LIBS)
c3vec.o:
	$(CC) src/C3Vec.cpp -c $(FLAGS) $(LIBS)
c3mat.o:
	$(CC) src/C3Mat.cpp -c $(FLAGS) $(LIBS)

clear:
	rm -f *o
mr_proper:
	rm -f *o $(EXEC)
