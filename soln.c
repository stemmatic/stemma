/*
	soln.c - driver for stemmatic inference
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "cache.h"
#include "hyp.h"
#include "boot.h"
#include "log.h"
#include "netmsg.h"

static Taxa *	readTaxa(Log *lg);

Length RetCost = 0;
Length MaxMP2 = 0;

int
	main(int argc, char *argv[])
{
	Taxa *taxa;
	Net *nt;
	Log *lg;
	FILE *fcon;
	Cache *curr;
	int present;
	char *cname, *ms;
	Length got0, got1;

	lg = logInit(argc, argv, 4, "[base [SOLN]]");

	taxa = readTaxa(lg);
	nt = ntNew(taxa);
	MaxMP2 = taxa->nVunits / 20 + 1;

	// Read constraints
	if ((fcon = logFile(lg, lgNO)) && ntConstraints(nt, fcon) != YES)
		return -1;
	fclose(fcon);

	curr = cacheNew(nt);
	cname = (argc == 2) ? "SOLN" : argv[2];
	present = cacheOpen(curr, cname);
	if (!present) {
		fprintf(stderr, "%s: no such cache.\n", cname);
		exit(-1);
	}
{ char *ret; if ((ret = getenv("RET"))) RetCost = atoi(ret); }
	cacheRead(nt, curr);
	cacheClose(curr);
	got0 = cacheRestore(nt, curr);
	got1 = ntCost(nt);

	if (got0 != got1 || getenv("FIX")) {
		ntMsg(nt, "Cached cost (%C) != current cost (%C), resetting cache...\n",
			got0, got1);
		cacheReset(curr);
		cacheSave(nt, got1, curr);
		cacheOpen(curr, cname);
		cacheLock(curr);
		cacheWrite(nt, curr);
		cacheUnlock(curr);
		cacheClose(curr);
		exit(-1);
	}
	
	if ((ms = getenv("MS")))
		logUncollate(lg, nt, ms);
	else if (!getenv("NOREPORT"))
		logResults(lg, nt);

	logAnalysis(stdout, nt);

	cacheFree(nt, curr);
	txFree(&taxa);

	exit(0);
	return 0;
}

// Read Taxa from file
static Taxa *
	readTaxa(Log *lg)
{
	FILE *fTaxa;
	Taxa *tx;

	fTaxa = logFile(lg, lgTX);
	if (!fTaxa) {
		fprintf(stderr, "Cannot open taxon-file\n");
		exit(-1);
		return (Taxa *) 0;
	}

	tx = txScan(fTaxa);
	fclose(fTaxa);
	if (!tx) {
		fprintf(stderr, "Fatal error reading taxon-file\n");
		exit(-1);
		return (Taxa *) 0;
	}

	return tx;
}
