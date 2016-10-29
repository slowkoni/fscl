
OBJ = fscl.o logmsg.o ms-input.o snp-input.o background-fsp.o sm-spline.o \
	sm-search.o scan-chromosome.o asc-bias.o cmdline-utils.o ms-parser.o ms-scanner.o
CC=gcc
CFLAGS=-Wall -ggdb -I ~/libs/include -m64 -O2 -march=native

LD_FLAGS=-m64 -L${HOME}/libs/lib -lm -lpthread -lgsl -lgslcblas

YACC=bison
LEX=flex

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

ms-parser.c: ms-parser.y
	$(YACC) -p ms -d -o$@ $<

ms-scanner.c: ms-scanner.lex ms-parser.c 
	$(LEX) -Cfa -Pms -o$@ $<

sm-sample: sm-spline.o sm-sample.o logmsg.o background-fsp.o asc-bias.o cmdline-utils.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

ascbias-segments: ascbias-segments.o cmdline-utils.o ms-parser.o ms-scanner.o
	$(CC) -o $@ $^ $(LD_FLAGS) $(LIBS) -lgsl -lgslcblas

clean:
	rm -f *.o *~ fscl ascbias-segments sm-sample
