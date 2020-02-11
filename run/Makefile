include ../*.mk

## check and recompile all *.o files upon which ALL depends
$(ODIR)/%.o: $(SDIR)/%.cpp $(ALLDEPS)
	$(CC) -c -o $@ $< $(AFLAGS) $(CFLAGS)

.PHONY: all
all: $(MAIN) $(PARAMETERS)

## link *.o files to create MAIN
$(MAIN): $(MAINOBJ)
	$(CC) -o $@ $^ $(AFLAGS) $(LIBS)

## link *.o files to create MAIN
$(PARAMETERS): $(PARAMETERSOBJ)
	$(CC) -o $@ $^ $(AFLAGS) $(LIBS)

## prevent creation of CLEAN file
.PHONY: cleanall cleanobj clean~ cleanprog

## clean directories of *.o and ~ files
cleanall: cleanobj clean~ cleanprog

cleanobj:
	-rm $(ODIR)/*.o

clean~:
	-rm $(SDIR)/*~ $(IDIR)/*~

cleanprog:
	-rm $(MAIN) $(PARAMETERS)