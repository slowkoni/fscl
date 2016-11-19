
OBJ = fscl.o logmsg.o ms-input.o snp-input.o background-fsp.o sm-spline.o \
	sm-search.o scan-chromosome.o asc-bias.o cmdline-utils.o ms-parser.o \
	ms-scanner.o pjh.o

CC=gcc
CFLAGS?=-Wall -ggdb -I . -m64 -O2 -march=native -fopenmp -DCUDA

LD_FLAGS?=-m64 -lm -lpthread -lgsl -lgslcblas -lcudart -lgomp

NVCC=/usr/local/cuda/bin/nvcc
NVCC_FLAGS?=-O2 -g

YACC=bison
LEX=flex

OS=$(shell uname)
ifeq ($(OS),Linux)
	CFLAGS+=-DLINUX
endif
ifeq ($(OS),Darwin)
	CFLAGS+=-DDARWIN
endif

all: fscl #sm-sample ascbias-segments

%.o: %.c
	$(CC) $(CFLAGS) -c $<

%.o: %.cu
	$(NVCC) $(NVCC_FLAGS) -c $<

fscl: $(OBJ)
	gcc -o $@ $^ $(LD_FLAGS) $(LIBS)

ms-parser.c: ms-parser.y
	$(YACC) -p ms -d -o$@ $<

ms-scanner.o: ms-scanner.c
	$(CC) $(CFLAGS) -Wno-unused-function -c $<

ms-scanner.c: ms-scanner.lex ms-parser.c 
	$(LEX) -Cfa -Pms -o$@ $<

sm-sample: sm-spline.o sm-sample.o logmsg.o background-fsp.o asc-bias.o cmdline-utils.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

ascbias-segments: ascbias-segments.o cmdline-utils.o ms-parser.o ms-scanner.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

clean:
	rm -f *.o *~ fscl ascbias-segments sm-sample
