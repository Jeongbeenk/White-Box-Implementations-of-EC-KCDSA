# Makefile - WBECDSA 실행을 위한 makefile

CFLAGS   = -std=c17
CFLAGS  += -I/opt/homebrew/Cellar/gmp/6.3.0/include
LDFLAGS  = -L/opt/homebrew/Cellar/gmp/6.3.0/lib
LDFLAGS += -lgmp -lm

SOURCES = main.c ge.c params.c utils.c EC.c sha256.c

.PHONY: all clean

all: main

main: $(SOURCES)
	$(CC) -o main $(SOURCES) $(CFLAGS) $(LDFLAGS)

clean:
	rm -f main
