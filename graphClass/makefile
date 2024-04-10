CC=/usr/bin/clang++

FLAGS = -g3 -DLINUX -fmessage-length=0 -std=c++11 -Wno-deprecated-register 
#FLAGS =  -O1 -g -fsanitize=address -fno-omit-frame-pointer  -DLINUX -fmessage-length=0 -std=c++11
#FLAGS = -g3 -fsanitize=address -fno-omit-frame-pointer  -DLINUX -fmessage-length=0 -std=c++11 
 
.cpp.o:
	$(CC) -c $(FLAGS) $<

.cc.o:
	$(CC) -c $(FLAGS) $<

.c.o:
	/usr/bin/clang -c $(FLAGS) $<

OBJS= main.o

all: Graph

Graph:$(OBJS)
	/usr/bin/clang++ *.o -L -lLAPACK  -lpthread -lm -o Graph
#	/usr/bin/clang++ *.o -static-libgcc -g -fsanitize=address -L -lLAPACK -lpthread -lm -o Graph
#	/usr/bin/clang++ *.o -static-libgcc -g -fsanitize=address -L -lLAPACK -lpthread -lm -o Graph
#	/usr/bin/clang++ *.o -static-libgcc -O -L -lLAPACK -lpthread -lm -lmpich -o Graph

