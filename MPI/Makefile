CXX=mpiicpc
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -pedantic -std=c++11 -O3 # -Werror
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=main.o grid.o double_buffer.o simulation.o

all: MPI 

MPI: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -rf hello MPI *.o *~
