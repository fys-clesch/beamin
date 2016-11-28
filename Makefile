CC = gcc

ifndef DEBUG
DEBUG = 0
endif

# 1 is the identifier for windows, 0 for linux:
ifndef OSid
OSid = 1
endif

ifeq ($(DEBUG),1)
CFLAGS = -O0 -g -Wall -Wextra -pedantic -lm -std=c99
else ifeq ($(DEBUG),2)
CFLAGS = -O0 -g -Wall -Wextra -pedantic -lm -std=c99 -fprofile-arcs -ftest-coverage
else
CFLAGS = -O3 -lm -std=c99
endif

OBJECTS = Main.o allocer.o nonlin_fit.o optic.o msg.o ff.o

GPF = Main.gcno allocer.gcno nonlin_fit.gcno optic.gcno msg.gcno ff.gcno

MIN: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS)
ifeq ($(OSid),1)
	@echo.
else
	@echo
endif
	@echo make done

# build object files:
%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
ifeq ($(OSid),1)
	del $(OBJECTS)
ifeq ($(DEBUG),2)
	del $(GPF)
endif
	@echo.
	@echo cleaned
else
	rm $(OBJECTS)
ifeq ($(DEBUG),2)
	rm $(GPF)
endif
	@echo
	@echo cleaned
endif
