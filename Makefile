bin/linkphase:	src/*
	g++ -o bin/linkphase src/linkphase.cpp -lboost_program_options

test: bin/linkphase
	./test.sh

clean:
	rm bin/*