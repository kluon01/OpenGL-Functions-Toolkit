# define g++ flags
CC = g++ -Wall -Wno-deprecated-declarations -O3

project4: project4.cpp
	$(CC) -o project4 project4.cpp $(LIB)
