Library			= 	libgswteos-10
Program			= 	gsw_check

PROGRAM_SOURCES = 	gsw_check_functions.c\
				  	gsw_oceanographic_toolbox.c\
				  	gsw_saar.c
LIBRARY_SRCS	= 	gsw_oceanographic_toolbox.c \
					gsw_saar.c

all: $(Program) $(Library)

$(Program):
	$(CPP) $(CRELEASE) $(PROGRAM_SOURCES) $(OUT)$(Program)$(X)

$(Library):
	$(CPP) $(LFLAGS) $(LIBRARY_SRCS) $(OUT)$(Library)	

clean:
    $(RM) *.exe *.o *.obj *.ilk *.pdb *.tmp *.i *~
	$(RM) $(Library)$(SHARED_POSTFIX)