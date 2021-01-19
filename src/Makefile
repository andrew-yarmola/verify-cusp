CC := g++ -Wall
CXX := g++ -Wall
COMMON := ACJ.o box.o elimination.o SL2ACJ.o Complex.o roundoff.o

verify: verify.o $(COMMON)
	$(CC) verify.o $(COMMON) -lm -o ../bin/verify

identify: identify.o $(COMMON)
	$(CC) identify.o $(COMMON) -lm -o ../bin/identify

rootcat: rootcat.o
	$(CC) rootcat.o -o ../bin/rootcat

# never optimize here, since the timing of arithmetic is important
test_float: test_float.o roundoff.o roundoff.h
	$(CC) test_float.o roundoff.o -o ../bin/test_float

