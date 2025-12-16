CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c17 -O3 -march=native -flto
CFLAGS_DEBUG = -Wall -Wextra -std=c17 -g -O0 -fsanitize=address,undefined
LDFLAGS = -lm

SRC = $(wildcard src/**/*.c src/*.c)
OBJ = $(SRC:.c=.o)
TARGET = debris

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

debug: CFLAGS = $(CFLAGS_DEBUG)
debug: $(TARGET)_debug

$(TARGET)_debug: $(OBJ)
	$(CC) $(CFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(OBJ) $(TARGET) $(TARGET)_debug

test:
	# Your test runner will go here

.PHONY: all clean debug test