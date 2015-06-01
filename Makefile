ifndef CC
CC=gcc
endif

CFLAGS:= $(CFLAGS) -Wall
PERF=-O2
DBG=-DDEBUG -O0
OBJDIR=./out
SRCDIR=./src
BINDIR=./bin
EXECUTABLES=simphy simphy_dbg simphy_static simphy_sorthologs

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

$(BINDIR)/simphy: $(SRCDIR)/main.c $(OBJECTS)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nSimphy built"

$(BINDIR)/simphy_sorthologs: $(SRCDIR)/main.c $(OBJECTS)
	mkdir -p $(BINDIR)
	@echo "\nThis target has been designed for internal usage and may not work properly in your system\n"
	$(CC) -D SORTHOLOGS $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nSimphy_sorthologs built"

$(BINDIR)/simphy_dbg: $(SRCDIR)/main.c $(DBG_OBJECTS)
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) -g $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nDebug version of Simphy built"

$(BINDIR)/simphy_static: $(SRCDIR)/main.c $(OBJECTS)
	mkdir -p $(BINDIR)
	@echo "\nThis target has been designed for internal usage and may not work properly in your system\n"
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(LIBS)
	@echo "\nStatic-linked version of Simphy built"

$(OBJDIR)/num_methods.o: $(SRCDIR)/num_methods.c $(SRCDIR)/num_methods.h $(SRCDIR)/trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/num_methods_dbg.o: $(SRCDIR)/num_methods.c $(SRCDIR)/num_methods.h $(SRCDIR)/trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/sampling.o: $(SRCDIR)/sampling.c $(SRCDIR)/sampling.h $(SRCDIR)/trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/sampling_dbg.o: $(SRCDIR)/sampling.c $(SRCDIR)/sampling.h $(SRCDIR)/trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/sql_managing.o: $(SRCDIR)/sql_managing.c $(SRCDIR)/sql_managing.h $(SRCDIR)/trees.h $(SRCDIR)/sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/sql_managing_dbg.o: $(SRCDIR)/sql_managing.c $(SRCDIR)/sql_managing.h $(SRCDIR)/trees.h $(SRCDIR)/sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

$(OBJDIR)/trees.o: $(SRCDIR)/trees.c $(SRCDIR)/trees.h $(SRCDIR)/sampling.h $(SRCDIR)/num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $< -c -o $@

$(OBJDIR)/trees_dbg.o: $(SRCDIR)/trees.c $(SRCDIR)/trees.h $(SRCDIR)/sampling.h $(SRCDIR)/num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@

.PHONY: clean
.PHONY: all
.PHONY: debug
.PHONY: static
.PHONY: sorthologs
.PHONY: simphy
.PHONY: simphy_debug
.PHONY: simphy_static
.PHONY: simphy_sorthologs

simphy: $(BINDIR)/simphy
simphy_debug: $(BINDIR)/simphy_dbg
simphy_dbg: simphy_debug
simphy_static: $(BINDIR)/simphy_static
simphy_sorthologs: $(BINDIR)/simphy_sorthologs
sorthologs: simphy_sorthologs
static: simphy_static
debug: simphy_dbg
all: $(EXECUTABLES)

clean:
	@echo "Cleaning object files directory\n"
	rm -f $(OBJECTS) $(DBG_OBJECTS)
	rm -rf $(OBJDIR)
	rm -rf simphy_dbg.dSYM
