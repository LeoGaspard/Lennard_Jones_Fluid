CC = g++
EXEC = Project
LIBS = -I /usr/local/boost_1_74_0
FLAGS = -fopenmp -Wall -std=c++2a

all: main.o
	$(CC) *.o -o $(EXEC) $(LIBS) $(FLAGS)
main.o : cdynamic.o
	$(CC) main.cpp -c $(FLAGS) $(LIBS)
cdynamic.o : capplication.o
	$(CC) CDynamic.cpp -c $(FLAGS) $(LIBS)
capplication.o: cbox.o
	$(CC) CApplication.cpp -c $(FLAGS) $(LIBS)
cbox.o: catom.o c3vec.o
	$(CC) CBox.cpp -c $(FLAGS) $(LIBS)
catom.o: cpos.o cspeed.o cforce.o
	$(CC) CAtom.cpp -c $(FLAGS) $(LIBS)
cpos.o: c3vec.o
	$(CC) CPos.cpp -c $(FLAGS) $(LIBS)
cspeed.o: c3vec.o
	$(CC) CSpeed.cpp -c $(FLAGS) $(LIBS)
cforce.o: c3vec.o
	$(CC) CForce.cpp -c $(FLAGS) $(LIBS)
c3vec.o:
	$(CC) C3Vec.cpp -c $(FLAGS) $(LIBS)

clear:
	rm -f *o
mr_proper:
	rm -f *o $(EXEC)
