CC=gcc
CFLAGS:= $(CFLAGS) -Wall
DBG=-DDEBUG
OBJDIR=out


EXECUTABLES=simphy simphy_dbg
_OBJECTS= num_methods.o sampling.o sql_managing.o trees.o
OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_OBJECTS))
_DBG_OBJECTS= num_methods_dbg.o sampling_dbg.o sql_managing_dbg.o trees_dbg.o
DBG_OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_DBG_OBJECTS))
#LIBS=-lm -lgsl -lcblas -lsqlite3 Improves the performance but reduces the compatibility.
LIBS=-lm -lgsl -lgslcblas -lsqlite3 -lmpfr

simphy: main.c $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)
	@echo "\nSimphy built"

simphy_dbg: main.c $(DBG_OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) -g $^ -o $@ $(LIBS)
	@echo "\nDebug version of Simphy built"

$(OBJDIR)/num_methods.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -c -o $@ $(LIBS)

$(OBJDIR)/num_methods_dbg.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@ $(LIBS)

$(OBJDIR)/sampling.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -c -o $@ $(LIBS)

$(OBJDIR)/sampling_dbg.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@ $(LIBS)

$(OBJDIR)/sql_managing.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -c -o $@ $(LIBS)

$(OBJDIR)/sql_managing_dbg.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@ $(LIBS)

$(OBJDIR)/trees.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -c -o $@ $(LIBS)

$(OBJDIR)/trees_dbg.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) $< -c -o $@ $(LIBS)

.PHONY: clean

clean:
	@echo "Cleaning object files directory\n"
	rm -f $(OBJECTS) $(DBG_OBJECTS)
	rmdir $(OBJDIR)
	rm -rf simphy_dbg.dSYM
