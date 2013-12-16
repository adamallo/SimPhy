CC=gcc
CFLAGS=-Wall
DBG=-DDEBUG
OBJDIR=out


EXECUTABLES=simphy simphy_dbg
_OBJECTS= num_methods.o sampling.o sql_managing.o trees.o
OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_OBJECTS))
_DBG_OBJECTS= num_methods_dbg.o sampling_dbg.o sql_managing_dbg.o trees_dbg.o
DBG_OBJECTS=$(patsubst %,$(OBJDIR)/%,$(_DBG_OBJECTS))
LIBS=-lm -lmpfr -lgsl -lcblas -lsqlite3

simphy: main.c $(OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@
	@echo "\nSimphy built"

simphy_dbg: main.c $(DBG_OBJECTS)
	$(CC) $(CFLAGS) $(DBG) $(LIBS) -g $^ -o $@
	@echo "\nDebug version of Simphy built"

$(OBJDIR)/num_methods.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $< -c -o $@

$(OBJDIR)/num_methods_dbg.o: num_methods.c num_methods.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $(DBG) $< -c -o $@

$(OBJDIR)/sampling.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $< -c -o $@

$(OBJDIR)/sampling_dbg.o: sampling.c sampling.h trees.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $(DBG) $< -c -o $@

$(OBJDIR)/sql_managing.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $< -c -o $@

$(OBJDIR)/sql_managing_dbg.o: sql_managing.c sql_managing.h trees.h sampling.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $(DBG) $< -c -o $@

$(OBJDIR)/trees.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $< -c -o $@

$(OBJDIR)/trees_dbg.o: trees.c trees.h sampling.h num_methods.h
	mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LIBS) $(DBG) $< -c -o $@


#test: indelap alignedPeps.phylip
#	./indelap -a alignedPeps.phylip
#	@echo
#	@echo "This test only verified correct running of program, not correct results."

.PHONY: clean

clean:
	@echo "Cleaning object files directory\n"
	rm -f $(OBJECTS) $(DBG_OBJECTS)
	rmdir $(OBJDIR)
	rm -rf simphy_dbg.dSYM
