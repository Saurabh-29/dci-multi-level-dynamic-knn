
output: example.o utils.o dci.o
	g++ -std=c++11 example.o utils.o dci.o -o output  -L/usr/local/lib -m64 -flto -fopenmp -dciENMP -O3 -Iinclude -lm -L/usr/lib/libblas -Wl,-rpath /usr/lib/libblas -lblas

utils.o: utils.cpp utils.h
	g++ -c -std=c++11 utils.cpp  -L/usr/local/lib -m64 -flto -fopenmp -dciENMP -O3 -Iinclude -lm -L/usr/lib/libblas -Wl,-rpath /usr/lib/libblas -lblas

dci.o: dci.cpp dci.h
	g++ -c -std=c++11 dci.cpp  -L/usr/local/lib -m64 -flto -fopenmp -dciENMP -O3 -Iinclude -lm -L/usr/lib/libblas -Wl,-rpath /usr/lib/libblas -lblas

example.o: example.cpp
	g++ -c -std=c++11 example.cpp  -L/usr/local/lib -m64 -flto -fopenmp -dciENMP -O3 -Iinclude -lm -L/usr/lib/libblas -Wl,-rpath /usr/lib/libblas -lblas

clean:
	rm *.o output