
OBJ = fscl.o logmsg.o ms-input.o snp-input.o background-fsp.o sm-spline.o \
	sm-search.o scan-chromosome.o asc-bias.o
CC=gcc
CFLAGS=-Wall -ggdb -I ~/libs/include -m64 -O2 -march=native

LD_FLAGS=-m64 -L${HOME}/libs/lib -lcmdline -lmsparser -lpthread -lgsl -lgslcblas

OS=$(shell uname)
ifeq ($(OS),Linux)
	CFLAGS+=-DLINUX
endif
ifeq ($(OS),Darwin)
	CFLAGS+=-DDARWIN
endif


%.o: %.c
	$(CC) $(CFLAGS) -c $<

all: fscl sm-sample ascbias-segments

fscl: $(OBJ)
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS)

sm-sample: sm-spline.o sm-sample.o logmsg.o background-fsp.o asc-bias.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

ascbias-segments: ascbias-segments.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

clean:
	rm -f *.o *~ fscl
