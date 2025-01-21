
/*
	cache.c - routines for caching/checkpointing solutions
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <sys/file.h>
#include <errno.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "netmsg.h"
#include "cache.h"

// Opaque type
struct Cache {
	FILE *fp;
	enum verbosity verbose;	// Verbose level
	Length cost;			// Cost of this solution
	Length obsolete;		// Time debt of this solution [OBSOLETE]
	Length rootCost;		// Cost of root node
	Flag cached;			// Cached?

	int nLinks;				// Current number of links
	Link *work, *end;		// Store list of links

	Cursor maxTax;			// Max Taxon in use + 1
	Flag *inuses;			// Which are in use.
	int nTaxa;				// Number of taxa in use
	int nNodes;				// Number of nodes in use == nTaxa 
	Cursor outgroup;		// Outgroup taxon.

	unsigned long *codes;	// Serial numbers of taxa
	Stratum *time;			// Time strata [OBSOLETE]
	vunit **states;			// States of taxa
#if DO_POLE
	Length *poles;			// [T] Precomputed Pole Positions
#endif
	Flag **noanc;			// Constraint of no ancestors [OBSOLETE]
};

////////////////////////////////////////////////////////

/* Quickie helpers
*/

int
	cacheLinks(Cache *cache)
{
	return cache->nLinks;
}

int
	cacheNodes(Cache *cache)
{
	return cache->nTaxa;
}

Length
	cacheCost(Cache *cache)
{
	return cache->cost;
}

Length
	cacheRootCost(Cache *cache)
{
	return cache->rootCost;
}

void
	cacheReset(Cache *cache)
{
	cache->cost = -1;
	cache->obsolete = -1;
	cache->rootCost = -1;
	cache->cached = NO;
}

void
	cacheMark(Cache *cache)
{
	cache->cached = NO;
}

int
	cacheCached(Cache *cache)
{
	return cache->cached;
}

enum verbosity
	cacheVerbosity(Cache *cache)
{
	return cache->verbose;
}

Link *
	cacheListLinks(Cache *cache, Link **work)
{
	*work = cache->work;
	return cache->end;
}

int
	cacheMixedNodes(Cache *cache)
{
	int *nParents, ii, nMixed=0;
	Link *link;

	newmem(nParents, cache->maxTax);
	ASET(nParents, 0, cache->maxTax);

	for (link = cache->work; link < cache->end; link++)
		nParents[link->to]++;
	for (ii = 0; ii < cache->maxTax; ii++)
		if (nParents[ii] > 1) nMixed++;

	free(nParents);
	return nMixed;
}


//////////////////////////////////
//
//	Cache and Restore: Saves solution states
//
Cache *
	cacheNew(Net *nt)
{
	register Taxa *tx = nt->taxa;
	Cache *cache;
	Cursor tt;

	newmem(cache, 1);

	cache->fp = 0;
	cache->verbose = C_VPASS;
	cache->nLinks = 0;
	cache->nNodes = 0;
	cache->end = cache->work = (Link *) 0;
	cacheReset(cache);

	newmem(cache->codes, tx->nTotal);
	newmem(cache->time, tx->nTotal);
	newmem(cache->states, tx->nTotal);
#if DO_POLE
	newmem(cache->poles, tx->nTotal);
	ZERO(cache->poles, tx->nTotal);
#endif

	NEWMAT(cache->noanc, tx->nTotal, tx->nTotal);
	ASET(cache->noanc[0], 0, tx->nTotal * tx->nTotal);

	cache->maxTax = nt->maxTax;
	newmem(cache->inuses, tx->nTotal);
	cache->nTaxa = nt->nTaxa;
	cache->outgroup = -1;

	ASET(cache->codes,~0,tx->nTotal);
	ZERO(cache->time,tx->nTotal);
	ZERO(cache->inuses,tx->nTotal);

	for (tt = 0; tt < tx->nTotal; tt++) {
		newmem(cache->states[tt], tx->nVunits);
		ASET(cache->states[tt], MISSING, tx->nVunits);
	}

	return cache;
}

void
	cacheFree(Net *nt, Cache *cache)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;

	cacheClose(cache);

	FREMAT(cache->noanc);
	for (tt = 0; tt < tx->nTotal; tt++) {
		free(cache->states[tt]);
	}

#if DO_POLE
	free(cache->poles);
#endif
	free(cache->inuses);
	free(cache->states);
	free(cache->time);
	free(cache->codes);
	free(cache->work);
	free(cache);
}

int
	cacheBetter(Net *nt, Length cost, Length rootCost, int nLinks,
		Cache *cache)
{
	register Taxa *tx = nt->taxa;

	if (cost > cache->cost)
		return NO;
	if (cost < cache->cost)
		return YES;

	if (rootCost == -1)
		rootCost = txApographic(tx,ntRoot(nt),nt->outgroup);
	if (rootCost < cache->rootCost)
		return YES;
	if (rootCost > cache->rootCost)
		return NO;

#if 1
	if (nLinks < cache->nLinks)
		return YES;
#endif

	return NO;
}

int	
	cacheComp(Net *nt, Cache *a, Cache *b)
{
	return cacheBetter(nt, a->cost, a->rootCost, a->nLinks, b);
}

#if DO_POLE
int
	cachePoleConstrained(Net *nt)
{
	register Taxa *tx = nt->taxa;
	Cursor to;

	for (to = 0; to < nt->maxTax; to++) {
		if (!nt->inuse[to])
			continue;
		if (nt->nParents[to] < 2)
			continue;

#if UNMIX
		for (Cursor *un = nt->unMixed; un < nt->endMixed; un++) {
			if (!nt->descendents[to][*un])
				continue;
			printf("\nMix ban violation for %d:'%s' --> %d:'%s'\n",
				to, txName(tx,to), *un, txName(tx,*un));
			fflush(stdout);
			return YES;
		}
#endif
		if (nt->banMixed[to]) {
			printf("\nMix ban violation for %d:'%s'\n", to, txName(tx,to));
			fflush(stdout);
			return YES;
		}
		if (!ntPoleCheck(nt, to)) {
			printf("\nPole problem for %d:'%s'\n", to, txName(tx,to));
			fflush(stdout);
			return YES;
		}
	}
	return NO;
}
#endif

int
	cacheSave(Net *nt, Length cost, Cache *cache)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;
	Length rootCost;

	// Quick check to avoid slower ntTime() and txApographic()
	if (cost > cache->cost)
		return NO;

	rootCost = txApographic(tx,ntRoot(nt),nt->outgroup);
	if (!cacheBetter(nt, cost, rootCost, nt->nLinks, cache))
		return NO;

#if 0 && DO_POLE
	// Very slow, so hoist the check closer to where it happens
	// until it never trips, then turn it off.
	assert( !cachePoleConstrained(nt) );
#endif

	cache->cached = YES;

	cache->cost = cost;
	cache->rootCost = rootCost;

	// Save connects in an array.
	if (cache->nLinks != nt->nLinks) {
		free(cache->work);
		cache->work = (Link *) 0;
	}
	cache->nLinks = nt->nLinks;
	cache->end = ntPreorderLinks(nt, &cache->work);

	// Save state
	cache->maxTax = nt->maxTax;
	cache->nNodes = 0;
	cache->nTaxa = nt->nTaxa;
	cache->outgroup = nt->outgroup;

	ACPY(cache->inuses,nt->inuse,tx->nTotal);
#if DO_POLE
	ACPY(cache->poles,nt->poles,tx->nTotal);
#endif

	for (tt = 0; tt < tx->nTotal; tt++) {
		if (nt->codes[tt] == cache->codes[tt])
			continue;
		cache->codes[tt] = nt->codes[tt];
		ACPY(cache->states[tt], txBase(tx,tt), tx->nVunits);
	}

	return YES;
}

void
	cacheNodeRestore(Net *nt, Cursor tt, Cache *cache)
{
	register Taxa *tx = nt->taxa;

	nt->inuse[tt] = cache->inuses[tt];

	if (nt->codes[tt] != cache->codes[tt]) {
		nt->codes[tt] = cache->codes[tt];
		ACPY(txBase(tx,tt), cache->states[tt], tx->nVunits);
#if DO_POLE
		if (tt >= tx->nExtant)
			nt->poles[tt] = cache->poles[tt];
#endif
		txPermuteTaxon(tx,tt);
	}
}

Length
	cacheRestore(Net *nt, Cache *cache)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;

	// Restore State
	nt->maxTax = cache->maxTax;
	nt->nTaxa = cache->nTaxa;
	nt->outgroup = cache->outgroup;

	ACPY(nt->inuse,cache->inuses,tx->nTotal);
#if DO_POLE
	ACPY(nt->poles,cache->poles,tx->nTotal);
#endif

	for (tt = 0; tt < tx->nTotal; tt++) {
		if (nt->codes[tt] != cache->codes[tt]) {
			nt->codes[tt] = cache->codes[tt];
			if (tt >= tx->nExtant) {
				ACPY(txBase(tx,tt), cache->states[tt], tx->nVunits);
				txPermuteTaxon(tx,tt);
			}
		}
	}
	// Restore connections
	ntRestoreLinks(nt, cache->work, cache->end);
	return cache->cost;
}

int
	cacheOpen(Cache *cache, char *cacheFile)
{
	FILE *fp;
	static char cacheBuf[128];
	int present = YES;

	strcpy(cacheBuf, (cacheFile || (cacheFile = getenv("CACHE")))
		? cacheFile : "BEST.cache");
	fp = fopen(cacheBuf, "r+b");
	if (!fp) {
		if (errno == ENOENT)
			fp = fopen(cacheBuf, "w+b");
		if (!fp) {
			fprintf(stderr, "Cannot open %s for caching.\n", cacheBuf);
			perror(cacheBuf);
			abort();
		}
		present = NO;
	}
	cache->fp = fp;
	return present;
}

int
	cacheLock(Cache *cache)
{
	assert( cache->fp != 0 );
	return flock(fileno(cache->fp), LOCK_EX) == OK;
}

int cacheUnlock(Cache *cache)
{
	assert( cache->fp != 0 );
	return flock(fileno(cache->fp), LOCK_UN) == OK;
}

int cacheClose(Cache *cache)
{
	if (cache->fp) {
		cacheUnlock(cache);
		fclose(cache->fp);
	}
	cache->fp = 0;
	return YES;
}

#define FREAD_CK(vp, sz, nn, fp) if (fread((vp), (sz), (nn), (fp)) != (nn)) abort(); else

// Return 1 if a better cached soln is read, 0 if a tie (and not read), -1 if the cached is worse.

int
	cacheRead(Net *nt, Cache *cache)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;
	int nTotal, nVunits, nLinks;
	unsigned int magic;
	Length cost, obsolete, rootCost, retCost;
	FILE *fp = cache->fp;
	size_t nb;

	assert( fp != 0 );
	rewind(fp);
	
	// Read Sanity Checks
	if (fread(&magic, sizeof magic, 1, fp) == 0)
		return -1;
	assert( magic == MAGIC );
	FREAD_CK(&nTotal, sizeof nTotal, 1, fp);
	assert( nTotal == tx->nTotal );
	FREAD_CK(&nVunits, sizeof nVunits, 1, fp);
	assert( nVunits == tx->nVunits );

	FREAD_CK(&nLinks, sizeof nLinks, 1, fp);

	FREAD_CK(&cost, sizeof cost, 1, fp);
	FREAD_CK(&obsolete, sizeof obsolete, 1, fp);
	FREAD_CK(&rootCost, sizeof rootCost, 1, fp);

	// Only Read if it is better.
	// Q:: call cacheBetter() instead?
	if (cost > cache->cost)
		return -1;
	if (cost == cache->cost) {
		if (rootCost > cache->rootCost)
			return -1;
		if (rootCost == cache->rootCost) {
			if (nt->nLinks > cache->nLinks)
				return -1;
			else if (nt->nLinks == cache->nLinks)
				return 0;
			// OK, we're good.
		}
	}

	cache->cost = cost;
	cache->rootCost = rootCost;

	// Read State
	FREAD_CK(&retCost, sizeof retCost, 1, fp);
if (!getenv("RET")) RetCost = retCost;
	FREAD_CK(&cache->maxTax, sizeof cache->maxTax, 1, fp);
	FREAD_CK(&cache->nTaxa, sizeof cache->nTaxa, 1, fp);
	FREAD_CK(&cache->nNodes, sizeof cache->nNodes, 1, fp);
	FREAD_CK(&cache->outgroup, sizeof cache->outgroup, 1, fp);
	FREAD_CK(cache->inuses, sizeof cache->inuses[0], tx->nTotal, fp);
	FREAD_CK(cache->time, sizeof cache->time[0], tx->nTotal, fp);
	for (tt = 0; tt < tx->nTotal; tt++) {
		cache->codes[tt] = TaxonCode++;
		FREAD_CK(cache->noanc[tt], sizeof cache->noanc[tt][0], tx->nTotal, fp);
		FREAD_CK(cache->states[tt], sizeof cache->states[tt][0], tx->nVunits, fp);
#if DO_POLE
		vunit *svBase = txBase(tx,tt);
		txBase(tx,tt) = cache->states[tt];
		cache->poles[tt] = txPole(tx, tt);
		txBase(tx,tt) = svBase;
#endif
	}

	// Read Links
	cache->nLinks = nLinks;
	if (cache->work)
		free(cache->work);
	newmem(cache->work, cache->nLinks);
	cache->end = cache->work + cache->nLinks;
	nb = fread(cache->work, sizeof cache->work[0], cache->nLinks, fp);

	assert( nb == cache->nLinks );

	return 1;
}

int
	cacheWrite(Net *nt, Cache *cache)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;
	unsigned int magic;
	FILE *fp = cache->fp;

	assert( fp != 0 );
	rewind(fp);
	
	// Write Sanity Checks
	magic = MAGIC;
	fwrite(&magic, sizeof magic, 1, fp);
	fwrite(&tx->nTotal, sizeof tx->nTotal, 1, fp);
	fwrite(&tx->nVunits, sizeof tx->nVunits, 1, fp);

	fwrite(&cache->nLinks, sizeof cache->nLinks, 1, fp);

	fwrite(&cache->cost, sizeof cache->cost, 1, fp);
	fwrite(&cache->obsolete, sizeof cache->obsolete, 1, fp);
	fwrite(&cache->rootCost, sizeof cache->rootCost, 1, fp);

	// Write State
	fwrite(&RetCost, sizeof RetCost, 1, fp);
	fwrite(&cache->maxTax, sizeof cache->maxTax, 1, fp);
	fwrite(&cache->nTaxa, sizeof cache->nTaxa, 1, fp);
	fwrite(&cache->nNodes, sizeof cache->nNodes, 1, fp);
	fwrite(&cache->outgroup, sizeof cache->outgroup, 1, fp);
	fwrite(cache->inuses, sizeof cache->inuses[0], tx->nTotal, fp);
	fwrite(cache->time, sizeof cache->time[0], tx->nTotal, fp);
	for (tt = 0; tt < tx->nTotal; tt++) {
		fwrite(cache->noanc[tt], sizeof cache->noanc[tt][0], tx->nTotal, fp);
		fwrite(cache->states[tt], sizeof cache->states[tt][0], tx->nVunits, fp);
	}

	// Write Links
	assert( cache->nLinks == (cache->end - cache->work) );
	fwrite(cache->work, sizeof cache->work[0], cache->nLinks, fp);
	fflush(fp);
	return YES;
}

////////////////////////////////////////////////////////////////////

void
	cacheSetVerbosity(Cache *cache, enum verbosity vlvl)
{
	cache->verbose = vlvl;
}

void
	cacheMsg(Net *nt, Cache *cache, enum verbosity vlvl, char *fmt, ...)
{
	va_list ap;

	// Don't call cacheMsg() with the upper and lower verbosity bounds
	assert( vlvl > C_VNONE && vlvl < C_VMAX );

	// Suppress messages if not wanted
	if (cache->verbose < vlvl)
		return;

	// Just call vprintf... and flush
	va_start(ap, fmt);
	ntVFMsg(stdout, nt, fmt, ap);
	va_end(ap);
	if (*fmt == '\r' || *fmt == '\n')
		fflush(stdout);
}

