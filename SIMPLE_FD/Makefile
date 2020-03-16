all: utils simple scratch

utils: utils.o
	g++ -c -o utils.o utils.cpp -O2

scratch: scratch_flow.o
	g++ -c -o scratch_flow.o scratch_flow.cpp -O2

simple: utils.o scratch_flow.o
	g++ -o scratch utils.o scratch_flow.o -O2
