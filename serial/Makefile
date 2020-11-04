CXX=icc
LD=${CXX}
CXXFLAGS+=-Wall -Wextra -pedantic -std=c++11 -O3 -Werror
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=main.o grid.o double_buffer.o simulation.o

all: serial

serial: $(OBJS)
	$(CXX)   $^ $(LDFLAGS+) -o $@

clean:
	rm -rf hello serial *.o *~
