Library			= 	libgswteos-10
Program			= 	gsw_check

PROGRAM_SOURCES = 	gsw_check_functions.c\
				  	gsw_oceanographic_toolbox.c\
				  	gsw_saar.c
LIBRARY_SRCS	= 	gsw_oceanographic_toolbox.c \
					gsw_saar.c

# \
!ifdef 0 #
include TOOLS.gcc # \
!endif # \
#

all: $(Program) $(Library)

$(Program):
		$(CPP) $(CRELEASE) $(PROGRAM_SOURCES) $(LIBS) $(OUT)$(Program)$(X)

library: $(Library)
$(Library):
		$(CPP) $(LFLAGS) $(LIBRARY_SRCS) $(LIBS) $(OUT)$(Library)$(SHARED_POSTFIX)

clean:
		$(RM) *.o *.obj *.ilk *.pdb *.tmp *.i *~
		$(RM) $(Library)$(SHARED_POSTFIX)
		$(RM) $(Program)$(X)
