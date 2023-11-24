CFLAGS = -std=c++11 -Wall -O3 -c -DDMEASURE #-DPROFILE_DEBUG #-DDEBUG #-DMICRO_DEBUG #-DPERTURB_CRIT #-DDEBUG_CHI_COMPUTATION
#CFLAGS = -std=c++11 -Wall -O0 -c -g -DDMEASURE -DCHI_DEBUG -DDEBUG #-DMICRO_DEBUG #-DPERTURB_CRIT #-DDEBUG_CHI_COMPUTATION

.PHONY = all

all: subgraph clean

subgraph: main.o query.o input.o vertex.o node.o
	g++ -o subgraph -Wall -O3 main.o query.o input.o vertex.o node.o

node.o: node.cpp node.h const.h
	g++ $(CFLAGS) node.cpp

vertex.o: vertex.cpp vertex.h node.h const.h
	g++ $(CFLAGS) vertex.cpp

query.o: query.cpp query.h node.h const.h
	g++ $(CFLAGS) query.cpp

input.o: input.cpp input.h vertex.h query.h
	g++ $(CFLAGS) input.cpp

main.o: main.cpp query.h input.h const.h
	g++ $(CFLAGS) main.cpp


clean:
	rm *.o
