#
#  $Id: Makefile,v 077b19612310 2013/09/17 22:55:03 fdelahoyde $
#
               CFLAGS:=	-O
            CINCLUDES:=	-I.
              Library:=	libgswteos-10.so
              Program:=	gsw_check_functions
      $(Program)_SRCS:=	gsw_check_functions.c \
			gsw_oceanographic_toolbox.c \
			gsw_saar.c
      $(Library)_SRCS:=	gsw_oceanographic_toolbox.c \
			gsw_saar.c
      $(Library)_OBJS:=	gsw_oceanographic_toolbox.o \
			gsw_saar.o

all:	$(Program)

$(Program):	$($(Program)_SRCS)
	gcc $(CFLAGS) $(CINCLUDES) -o $(Program) $($(Program)_SRCS) -lm

library:	$(Library)

$(Library):	$($(Library)_SRCS)
	gcc -fPIC -c $(CFLAGS) $(CINCLUDES) $($(Library)_SRCS)
	gcc -shared -o $(Library) $($(Library)_OBJS) -lm

clean:
	rm -f $(Program) $(Library) $($(Library)_OBJS)

