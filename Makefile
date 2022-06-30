bin/linkphase: src/linkphase.cpp src/pedigree.cpp src/AnimalInfo.cpp
	g++ -o bin/linkphase src/linkphase.cpp src/AnimalInfo.cpp -lboost_program_options

all: bin/linkphase

clean:
	rm bin/*

test: bin/linkphase
	./test.sh

