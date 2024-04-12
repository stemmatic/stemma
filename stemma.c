/*
	stemma.c - driver for stemmatic inference
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "cache.h"
#include "hyp.h"
#include "heur.h"
#include "boot.h"
#include "log.h"

static Taxa *	readTaxa(Log *lg);
static void		initStemma(Net *nt, Cache *curr, Cache *best);
static int		bootStemma(Net *nt, Cache *best, Log *lg, int reps);

Length RetCost = 0;
Length MaxMP2 = 0;

unsigned int Seed = 0;

int
	main(int argc, char *argv[])
{
	Taxa *taxa;
	Net *nt;
	Cache *curr, *best;
	Log *lg;
	FILE *fcon;
	int started, nX;
	char *verbose, *mode, *maxmix = 0;

	block {
		char *seed = getenv("SEED");
		if (!seed) {
			Seed = time((time_t *) 0);
			if ((seed = getenv("TAG")))
				Seed -= 100*atoi(seed);
		} else
			Seed = atoi(seed);
		srandom(Seed);
	}

	lg = logInit(argc, argv, 6, "{a(nneal),b(oot),r(atchet)");

	taxa = readTaxa(lg);
	nt = ntNew(taxa);
	MaxMP2 = taxa->nVunits / 20 + 1;

#if DO_MAXMIX
	if ((maxmix = getenv("MAXMIX")))
		nt->maxMix = atoi(maxmix);
#endif

	// Read constraints
	if ((fcon = logFile(lg, lgNO)) && ntConstraints(nt, fcon) != YES)
		return -1;
	fclose(fcon);

	curr = cacheNew(nt);
	best = cacheNew(nt);
	if ((verbose = getenv("VERBOSE"))) {
		cacheSetVerbosity(curr, atoi(verbose));
		cacheSetVerbosity(best, atoi(verbose));
	}

	cacheOpen(best, "START");
	cacheLock(best);
	started = cacheRead(nt, best);
	cacheUnlock(best);
	cacheClose(best);

	if (started >= 0) {
		cacheRestore(nt, best);
		cacheSave(nt, cacheCost(best), curr);
#if DO_MAXMIX
		if (RetCost == 0 && !maxmix)
			nt->maxMix = cacheMixedNodes(best);
#endif
	} else {
		if (maxmix) {
			RetCost = 0;
			initStemma(nt, curr, best);

			cacheRestore(nt, best);
			cacheSave(nt, cacheCost(best), curr);
		} else {
			printf("Determining reticulation cost...\n");

			RetCost = nt->taxa->nVunits;
			initStemma(nt, curr, best);

			cacheRestore(nt, best);
			cacheSave(nt, cacheCost(best), curr);

			// Set RetCost to average cost per node plus std dev.
			block {
				double mean = (double) cacheCost(best)/cacheNodes(best);
#if 0
				double dev, sumdev = 0.0;
				Cursor t;
				for (t = 0; t < nt->maxTax; t++) {
					if (!nt->inuse[t])
						continue;
					dev = nt->cumes[t] - mean;
					sumdev += dev * dev;
				}
				
				RetCost = ceil(mean + sqrt(sumdev/cacheNodes(best))/2.0);
#else
				RetCost = ceil(mean);
#endif
			}
			printf("Reticulation cost is: %d\n", (int) RetCost);
		}
		srandom(Seed);

		(void) ntCost(nt);
		cacheOpen(best, "START");
		cacheLock(best);
		cacheWrite(nt, best);
		cacheUnlock(best);
		cacheClose(best);
	}
{ char *ret; if ((ret = getenv("RET"))) RetCost = atoi(ret); }

	mode = (argv[2]) ? argv[2] : "";
	nX = nt->taxa->nExtant;
	if (*mode == 'b') {
		int nBoot = (argv[3]) ? atoi(argv[3]) : 1000;
		bootStemma(nt, best, lg, nBoot);
	} else if (*mode == 'r') {
		int nCycles = (argv[3]) ? atoi(argv[3]) : nX*nX*nX;
		int nRatchets = (argv[3] && argv[4]) ? atoi(argv[4]) : 1;
		heurRatchet(nt, curr, best, nCycles, nRatchets);
	} else if (*mode == 'a') {
		int temperature = (argv[3]) ? atoi(argv[3]) : nX*nX;
		heurAnneal(nt, curr, best, temperature);
	} else {
		/* Default to START sequence. */
		int startTemp = 2*nX;
		int startTimes = nX/2;
		int xxTemp = nX*nX;
		int xxTimes = nX*nX*nX;
		int i;

		switch (0) {
		default:
			if (!argv[2]) break;
			if (*argv[2] != '.') startTemp = atoi(argv[2]);
			if (!argv[3]) break;
			if (*argv[3] != '.') startTimes = atoi(argv[3]);
			if (!argv[4]) break;
			if (*argv[4] != '.') xxTemp = atoi(argv[4]);
			if (!argv[5]) break;
			if (*argv[5] != '.') xxTimes = atoi(argv[5]);
		}

		setenv("CACHE", "Sxx", 1);
		for (i = 0; i < startTimes; i++) {
			char tag[4];
			sprintf(tag, "%02d", i % 100);
			setenv("TAG", tag, 1);
			unlink("Sxx");
			cacheReset(best);
			cacheRestore(nt, curr);
			heurAnneal(nt, curr, best, startTemp);
		}
		cacheRestore(nt, best);
		unsetenv("TAG");
		heurAnneal(nt, curr, best, xxTemp);
		heurRatchet(nt, curr, best, xxTimes, 1);
	}

	cacheFree(nt, curr);
	cacheFree(nt, best);
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

// initStemma: generate initial stemma by greedy alg. and CRR

static void	
	initStemma(Net *nt, Cache *curr, Cache *best)
{
	Length cost;

	cost = ntCost(nt);
	cacheSave(nt, cost, curr);
	cacheSave(nt, cost, best);

	hypLinkUp(nt, curr, best);

	cost = hypCRR(nt, curr, best, "--> ");
}

static int	
	bootStemma(Net *nt, Cache *best, Log *lg, int reps)
{
	Boot *bt;

	cacheOpen(best, "SOLN");
	cacheLock(best);
	cacheRead(nt, best);
	cacheUnlock(best);
	cacheClose(best);

	bt = bootNew(nt, best);
	bootStrap(nt, bt, reps);
	logBoot(nt, bt, lg);
	bootFree(nt, bt);
	return YES;
}
