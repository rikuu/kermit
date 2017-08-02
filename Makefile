CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wextra -pedantic
CPPFLAGS=
INCLUDES=	-I. -Iminiasm
OBJS=		miniasm/sys.o miniasm/sdict.o miniasm/paf.o miniasm/asg.o miniasm/common.o miniasm/hit.o miniasm/asm.o cf.o color.o index.o marker.o
PROG=		kermit kermit-color
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

kermit:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

kermit-color:$(OBJS) color.o kermit-color.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

main.o: miniasm/kvec.h miniasm/sys.h miniasm/paf.h miniasm/sdict.h miniasm/miniasm.h miniasm/asg.h color.h
cf.o: miniasm/paf.h miniasm/kseq.h color.h
color.o: color.h miniasm/asg.h
