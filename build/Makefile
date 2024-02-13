FLAGS = -O3 -Wall -pedantic -std=c++17 -DNDBUG

lbm: main.o parameters.o lattice.o functions.o lbm.o obstacles.o
	g++ $(FLAGS) -o lbm main.o parameters.o lattice.o functions.o lbm.o obstacles.o

main.o: ../source/main.cpp ../include/parameters.h
	g++ $(FLAGS) -c ../source/main.cpp

parameters.o: ../source/parameters.cpp ../include/parameters.h
	g++ $(FLAGS) -c ../source/parameters.cpp

lattice.o: ../source/lattice.cpp ../include/lattice.h
	g++ $(FLAGS) -c ../source/lattice.cpp

functions.o: ../source/functions.cpp ../include/functions.h ../include/collections.h
	g++ $(FLAGS) -c ../source/functions.cpp

lbm.o: ../source/lbm.cpp ../include/lbm.h ../include/collections.h
	g++ $(FLAGS) -c ../source/lbm.cpp

obstacles.o: ../source/obstacles.cpp ../include/obstacles.h
	g++ $(FLAGS) -c ../source/obstacles.cpp

clean:
	rm -rf *.o lbm ../outputs/*.vtk
