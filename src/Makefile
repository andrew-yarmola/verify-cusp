CC=g++
LIBS=-lm
COMMON= ACJ.o box.o elimination.o params.o SL2ACJ.o Complex.o roundoff.o

.SUFFIXES: .o .c .C

verify: verify.o $(COMMON)
	$(CC) verify.o $(COMMON) $(LIBS) -o ../bin/verify

identify: identify.o $(COMMON)
	$(CC) identify.o $(COMMON) $(LIBS) -o ../bin/identify

rootcat: rootcat.o
	$(CC) rootcat.o -o ../bin/rootcat

# never optimize here, since the timing of arithmetic is important
test_float: test_float.o roundoff.o roundoff.h
	$(CC) test_float.o roundoff.o -o test_float

