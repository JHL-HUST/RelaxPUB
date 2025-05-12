CC=g++ -flto -O3
CFLAGS=-c -I. -std=c++11

all: RelaxGCB-kPlexS

RelaxGCB-kPlexS: .obj/main.o .obj/Graph.o .obj/Global.o
	${CC} .obj/main.o .obj/Graph.o .obj/Global.o -o RelaxGCB-kPlexS
	rm .obj/*.o

.obj/main.o: main.cpp
	${CC} ${CFLAGS} -o .obj/main.o main.cpp

.obj/Graph.o: Graph.cpp
	${CC} ${CFLAGS} -o .obj/Graph.o Graph.cpp

.obj/Global.o: Global.cpp
	${CC} ${CFLAGS} -o .obj/Global.o Global.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
