CFLAGS=-g -Wall 
#CFLAGS+=-O2       # enable optimization
CFLAGS+=-O3       # enable optimization
LDLIBS=-lm
MAKEFLAGS=--no-print-directory

# Enable profiling
#CFLAGS+=-pg
#LDFLAGS=-static -pg
#LDLIBS+=-lc_p

APPS=stemma soln

BIN=$(HOME)/bin

APPEXES=$(APPS:%=$(BIN)/%)
APPOBJS=taxon.o net.o netmsg.o hyp.o heur.o cache.o log.o boot.o state.o

OBJS=$(APPS:%=%.o) $(APPOBJS)

all:	$(APPEXES)
#	$(MAKE) -C Utest

$(APPEXES): $(BIN)/%: %
	ln -f $< $(BIN)/

$(APPS):	$(APPOBJS)

.PHONY:	ut
ut: $(APPOBJS)
	$(MAKE) -C Utest
	$(MAKE) -C Utest ut

.PHONY:	clean
clean:	
	-rm $(OBJS) $(APPS)
	$(MAKE) -C Utest clean

.PHONY:	cod
cod:
	rlog -R -L RCS/*

.PHONY:	cou
cou:
	co -u $(OBJS:.o=.c) $(OBJS:.o=.h)

.PHONY:	ci
ci:
	ci -u $(OBJS:.o=.c) $(OBJS:.o=.h)

.PHONY: depend
depend:
	makedepend -- $(CFLAGS) -- $(OBJS:.o=.c)

# DO NOT DELETE


boot.o: boot.h cache.h hyp.h net.h netmsg.h stemma.h taxon.h
cache.o: cache.h net.h netmsg.h stemma.h taxon.h
heur.o: cache.h heur.h hyp.h net.h netmsg.h state.h stemma.h taxon.h
hyp.o: boot.h cache.h hyp.h net.h netmsg.h state.h stemma.h taxon.h
log.o: boot.h cache.h log.h net.h netmsg.h stemma.h taxon.h
netmsg.o: net.h netmsg.h stemma.h taxon.h
net.o: net.h netmsg.h state.h stemma.h taxon.h
soln.o: boot.h cache.h hyp.h log.h net.h netmsg.h stemma.h taxon.h
state.o: net.h state.h stemma.h taxon.h
stemma.o: boot.h cache.h heur.h hyp.h log.h net.h stemma.h taxon.h
taxon.o: stemma.h taxon.h
