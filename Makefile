Library         = libgswteos-10
Program         = gsw_check

PROGRAM_SOURCES = gsw_check_functions.c\
                  gsw_oceanographic_toolbox.c\
                  gsw_saar.c
LIBRARY_SRCS    = gsw_oceanographic_toolbox.c \
                  gsw_saar.c 

# This includes TOOLS.gcc if make (unix) is used
# The #\ logic causes the include to be ignored by nmake (windows)
# nmake automatically includes TOOLS.ini
# \
!ifdef 0 #
include TOOLS.gcc # \
!endif # \
#

# --- NOTES ---:
# (1) The default CC variable is defined in the TOOLS.gcc file (make) / Tools.ini (nmake) file.
#     use (n)make CC=compiler_command_of_choice to override the default compiler.
# (2) On windows the default CC variable is cl with the option /TP. This causes the project to be compiled
#     as C++ code on windows. This is necessary because msvc does not support C99 style complex numbers.
#     See also notes in TOOLS.ini.


all: $(Program) $(Library)

$(Program):
		$(CC) $(CFLAGS) $(PROGRAM_SOURCES) $(LIBS) $(OUT)$(Program)$(X)

library: $(Library)
$(Library):
		$(CC) $(LFLAGS) $(LIBRARY_SRCS) $(LIBS) $(OUT)$(Library)$(SHARED_POSTFIX)

clean:
		$(RM) *.o *.obj *.ilk *.pdb *.tmp *.i *~
		$(RM) $(Library)$(SHARED_POSTFIX)
		$(RM) $(Program)$(X)

# Print the the CC variable (for debugging purpose).
# This is useful to check if the correct compiler is used with the specified commands.
name_compiler:
		@echo $(CC)

# Print the the EXTRA_FLAGS variable (for debugging purpose).
# This is useful to check if the correct extra flags are used with the specified commands.
name_extra_flags:
		@echo $(EXTRA_FLAGS)
