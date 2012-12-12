
LDFLAGS+=-lm -ggdb3
CFLAGS+= -Wextra -Wall -ggdb3


OBJECTS=igrf.o igrf-tst.o

all: igrf-tst tags

tags: $(wildcard *.c) $(wildcard *.h)
	ctags $(wildcard *.c) $(wildcard *.h)

igrf-tst: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -g -o $@

.PHONY:
clean:
	rm -vf igrf-tst
	rm -vf $(OBJECTS)



