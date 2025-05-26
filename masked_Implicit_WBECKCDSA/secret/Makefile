# Makefile - Random Encoding Matrix 생성을 위한 makefile
CC = gcc
CFLAGS   = -std=c11
CFLAGS  += -I/opt/homebrew/Cellar/gmp/6.3.0/include
LDFLAGS  = -L/opt/homebrew/Cellar/gmp/6.3.0/lib
LDFLAGS += -lgmp -lm 

SOURCES = main.c encode_gen.c ../public/EC.c secret_params.c ../public/params.c ../public/utils.c ../public/sha256.c

.PHONY: all clean

all: ki main

main: $(SOURCES)
	$(CC) -o main $(SOURCES) $(CFLAGS) $(LDFLAGS)

ki:
	# python3 k_gen.py

clean:
	rm -f main
