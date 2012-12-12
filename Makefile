
LDFLAGS+=-lm


OBJECTS=igrf.o igrf-tst.o

all: igrf-tst tags

tags: $(wildcard *.c) $(wildcard *.h)
	ctags $(wildcard *.c) $(wildcard *.h)

igrf-tst: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -g -o $@



