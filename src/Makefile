CC=gcc
CFLAGS:= $(CFLAGS) -Wall
PERF=-O2
DBG=-DDEBUG -O0
OBJDIR=out
EXECUTABLES=simphy simphy_dbg simphy_static

_OBJECTS= num_methods.o sampling.o sql_managing.o trees.o
OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_OBJECTS))
_DBG_OBJECTS= num_methods_dbg.o sampling_dbg.o sql_managing_dbg.o trees_dbg.o
DBG_OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_DBG_OBJECTS))

#-lcblas instead of lgslcblas Improves performance but reduces the compatibility.

#Dynamic libraries
C_LIBS= -lm -ldl -lpthread #Allways dynamically linked
D_LIBS= -lgsl -lgslcblas -lsqlite3 -lmpfr

#Static libraries for MAC
_S_LIBS= libgsl.a libgslcblas.a libsqlite3.a libgmp.a libmpfr.a
MS_LIBS=$(patsubst %,$(LD_LIBRARY_PATH)/%,$(_S_LIBS)) #BSD's LD needs the full path

#Static libraries for Linux
LS_LIBS= -lgsl -lgslcblas -lsqlite3 -lmpfr -lgmp

UNAME_S= $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LIBS= -Wl,-no_pie $(C_LIBS) $(MS_LIBS)
else
LIBS= -Wl,-Bstatic $(LS_LIBS) -Wl,-Bdynamic $(C_LIBS)
endif

simphy: main.c $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nSimphy built"

simphy_dbg: main.c $(DBG_OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) -g $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nDebug version of Simphy built"

simphy_static: main.c $(OBJECTS)
	@echo "\nThis target has been designed for internal usage and may not work properly in your system\n"
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(LIBS)
	@echo "\nStatic-linked version of Simphy built"

$(OBJDIR)/num_methods.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/num_methods_dbg.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/sampling.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/sampling_dbg.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/sql_managing.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/sql_managing_dbg.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/trees.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/trees_dbg.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

.PHONY: clean
.PHONY: all
.PHONY: debug
.PHONY: static

static: simphy_static
debug: simphy_dbg
all: $(EXECUTABLES)
clean:
	@echo "Cleaning object files directory\n"
	rm -f $(OBJECTS) $(DBG_OBJECTS)
	rmdir $(OBJDIR)
	rm -rf simphy_dbg.dSYM
