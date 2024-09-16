
/*
	hear.c - Heuristic searching routines
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdarg.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "netmsg.h"
#include "cache.h"
#include "hyp.h"
#include "heur.h"
#include "state.h"

#define FIRST_FOUND 1
#define DO_BRANCH 0
#define ABS_DELTA 1


/////////////////////////////////////////////////////////////
//
//	heurRatchet - Perturb by resampling, and Restore
//
//

Length
	heurRatchet(Net *nt, Cache *start, Cache *best,
		int nCycles, int nRatchets)
{
	Length got;
	int nC, nR;
	Cache *curr = cacheNew(nt);
	Cache *ratchet = cacheNew(nt);
	static char prefix[80];
	char *cache_file = getenv("CACHE");
	char *tag = getenv("TAG");
	int bump = getenv("BUMP") ? nCycles : 0;

	// Run the ratchet
	cacheSetVerbosity(curr, cacheVerbosity(best));
	cacheSetVerbosity(ratchet, cacheVerbosity(best));
	if (!tag)
		tag = "xx";
	if (!cache_file) {
		sprintf(prefix, "S%2s", tag);
		cache_file = prefix;
	}
	cacheOpen(best, cache_file);

	cacheSave(nt, cacheCost(start), curr);
	cacheSave(nt, cacheCost(start), best);
	cacheLock(best);
	if (cacheRead(nt, best) < 0)
		cacheWrite(nt, best);
	cacheUnlock(best);
	cacheRestore(nt, best);

	cacheMsg(nt, best, C_VPASS, "%2s ", tag);
	cacheMsg(nt, best, C_VPASS,
		"Ratchet: %6d of %6d / %2d ",
		0, nCycles, 1 - nRatchets);
	cacheMsg(nt, best, C_VPASS, "RET=%3d        ", RetCost);
	cacheMsg(nt, best, C_VPASS, "best: %C-%R ",
		cacheCost(best), cacheRootCost(best));
	cacheMsg(nt, best, C_VPASS, "u%.3f ",
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VPASS, "x%3d", cacheMixedNodes(best));
	cacheMsg(nt, best, C_VPASS, "\n");

	for (nC = 0; nC < nCycles; nC++) {
		for (nR = 0; nR < nRatchets; nR++) {
			// Phase 1: perturb and hill climb
			cacheReset(ratchet);
			txPerturb(nt->taxa, nt->codes);
			cacheReset(curr);
			got = ntCost(nt);
			cacheSave(nt, got, curr);
			cacheSave(nt, got, ratchet);
			sprintf(prefix, "*** [%d/%d, %d/%d: "CST_F"-"RTC_F"] *** ",
				nC+1, nCycles, nR+1, nRatchets,
				cacheCost(best), cacheRootCost(best));
			hypCRR(nt, curr, ratchet, prefix);
			cacheRestore(nt, ratchet);
			
			// Phase 2: restore and hill climb
			cacheReset(ratchet);
			txRestore(nt->taxa, nt->codes);
			got = ntCost(nt);
			cacheReset(curr);
			cacheSave(nt, got, curr);
			cacheSave(nt, got, ratchet);
			sprintf(prefix, "    [%d/%d, %d/%d: "CST_F"-"RTC_F"]     ",
				nC+1, nCycles, nR+1, nRatchets,
				cacheCost(best), cacheRootCost(best));
			hypCRR(nt, curr, ratchet, prefix);
			cacheRestore(nt, ratchet);

			cacheMark(best);
			cacheSave(nt, cacheCost(ratchet), best);

			cacheMsg(nt, best, C_VPASS, "%2s ", tag);
			cacheMsg(nt, best, C_VPASS,
				"Ratchet: %6d of %6d / %2d ",
				nC+1, nCycles, nR+1 - nRatchets);
			cacheMsg(nt, best, C_VPASS, "curr: %C-%R ",
				cacheCost(curr), cacheRootCost(curr));
			cacheMsg(nt, best, C_VPASS, "best: %C-%R ",
				cacheCost(best), cacheRootCost(best));
			cacheMsg(nt, best, C_VPASS, "u%.3f ",
				(double) cacheCost(best)/cacheNodes(best));
			cacheMsg(nt, best, C_VPASS, "x%3d", cacheMixedNodes(best));

			cacheMsg(nt, best, C_VPASS, "\r");
			cacheMsg(nt, best, (cacheCached(best)) ? C_VPASS : C_VITER, "\n");
			if (bump && cacheCached(best))
				nCycles = nC + bump;
		}
		cacheLock(best);
		if (cacheRead(nt, best) < 0)
			cacheWrite(nt, best);
		cacheUnlock(best);
		cacheRestore(nt, best);
		cacheSave(nt, cacheCost(best), best);
	}
	cacheMsg(nt, best, C_VRUN, "%2s best: %C-%R ",
		tag, cacheCost(best), cacheRootCost(best));
	cacheMsg(nt, best, C_VRUN, "u%.3f ",
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VRUN, "x%3d ", cacheMixedNodes(best));

	cacheOpen(best, "SOLN");
	cacheLock(best);
	if (cacheRead(nt, best) < 0)
		cacheWrite(nt, best);
	cacheUnlock(best);
	cacheClose(best);
	cacheMsg(nt, best, C_VRUN, "soln: %C-%R ",
		cacheCost(best), cacheRootCost(best));
	cacheMsg(nt, best, C_VRUN, "u%.3f ",
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VRUN, "x%3d", cacheMixedNodes(best));
	cacheMsg(nt, best, C_VRUN, "   \n");

	cacheFree(nt, curr);
	cacheFree(nt, ratchet);
	return cacheCost(best);
}

/////////////////////////////////////////////////////////////
//
//	heurAnneal - Simulated Annealing
//
//

Length
	heurAnneal(Net *nt, Cache *start, Cache *best, int initTemp)
{
	Length got;
	Cache *curr = cacheNew(nt);
	Cache *pool = cacheNew(nt);
	Cache *ratchet = cacheNew(nt);
	static char prefix[80];
	int temperature = initTemp;
	char *cache_file = getenv("CACHE");
	char *tag = getenv("TAG");
	char *pooling = getenv("POOL");
	time_t maxT=0;
	char *maxt;

	if ((maxt = getenv("MAXT")))
		maxT = time((time_t *) NULL) + atoi(maxt);

	cacheSetVerbosity(curr, cacheVerbosity(best));
	cacheSetVerbosity(pool, cacheVerbosity(best));
	cacheSetVerbosity(ratchet, cacheVerbosity(best));
	if (!tag)
		tag = "xx";
	if (!cache_file) {
		sprintf(prefix, "S%2s", tag);
		cache_file = prefix;
	}
	cacheOpen(best, cache_file);

	cacheSave(nt, cacheCost(start), curr);
	cacheSave(nt, cacheCost(start), best);
	cacheLock(best);
	if (cacheRead(nt, best) < 0)
		cacheWrite(nt, best);
	cacheUnlock(best);
	cacheRestore(nt, best);
	cacheSave(nt, cacheCost(best), curr);

	cacheMsg(nt, best, C_VPASS, "%2s ", tag);
	cacheMsg(nt, best, C_VPASS, "T%6d/%6d %C-%R u%.3f",
		temperature, 0, cacheCost(best), cacheRootCost(best),
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VPASS, " x%3d", cacheMixedNodes(best));
	cacheMsg(nt, best, C_VPASS, "         |: R=%ld              \n", RetCost);

	// Run the ratchet anneal
	do {
		int gotBest;
		int iter = 0;
		long delta;
		Length now, oldBest;
		double Boltz_prob = 0.0, rnd = 0.0;

#if !ABS_DELTA
		Length prev;
		prev = cacheCost(best);
#endif
		do {
			gotBest = NO;
			oldBest = cacheCost(best);
			iter++;

			// Phase 1: perturb and hill climb
			cacheReset(ratchet);
			txPerturb(nt->taxa, nt->codes);
			got = ntCost(nt);
			cacheReset(curr);
			cacheSave(nt, got, curr);
			assert( nt->nTaxa == cacheNodes(curr) );
			cacheSave(nt, got, ratchet);
			sprintf(prefix, "*** [T=%d, %d: "CST_F"-"RTC_F"] *** ",
				temperature, iter, cacheCost(best), cacheRootCost(best));
			hypCRR(nt, curr, ratchet, prefix);
			cacheRestore(nt, ratchet);
			
			// Phase 2: restore and hill climb
			cacheReset(ratchet);
			txRestore(nt->taxa, nt->codes);
			cacheReset(curr);
			got = ntCost(nt);
			cacheSave(nt, got, curr);
			cacheSave(nt, got, ratchet);
			sprintf(prefix, "    [T=%d, %d: "CST_F"-"RTC_F"]     ",
				temperature, iter, cacheCost(best), cacheRootCost(best));
			hypCRR(nt, curr, ratchet, prefix);
			cacheRestore(nt, ratchet);

			cacheSave(nt, cacheCost(ratchet), best);

			cacheMsg(nt, best, C_VITER, "\n");
			cacheMsg(nt, best, C_VPASS, "%2s ", tag);
			cacheMsg(nt, best, C_VPASS, "T%6d/%6d",
				temperature, iter);
			cacheMsg(nt, best, C_VPASS, " %C-%R ",
				cacheCost(curr), cacheRootCost(curr));
			now = cacheCost(ratchet);

			{
				Length currBest = cacheCost(best);
				if (currBest < oldBest) {
					delta = now - oldBest;
					oldBest = currBest;
					gotBest = YES;
				} else
					delta = now - currBest;
			}
#if !ABS_DELTA
			delta = now - prev;
			prev = now;
#endif
			if (maxt && time((time_t *)0) >= maxT) {
				delta = 100;
				temperature = 0.5;
			}
			if (delta >= 0) {
				// Transition based on Boltzman distribution
				if (delta > 0)
					Boltz_prob = exp(-(double)delta/temperature);
				else
					Boltz_prob = exp(-0.1/temperature);
				rnd = (double)random()/RAND_MAX;
				cacheMsg(nt, best, C_VPASS, "d%3ld, p%.3f >%.3f",
					delta, Boltz_prob, rnd);
				if (rnd >= Boltz_prob)
					cacheMsg(nt, best, C_VPASS, "!");
				else
					cacheMsg(nt, best, C_VPASS, " ");
			} else
				cacheMsg(nt, best, C_VPASS, "                    ");
			cacheMsg(nt, best, C_VPASS, " | %C-%R ",
				cacheCost(best), cacheRootCost(best));
			cacheMsg(nt, best, C_VPASS, "u%.3f ",
				(double) cacheCost(best)/cacheNodes(best));
			cacheMsg(nt, best, C_VPASS, "x%3d", cacheMixedNodes(best));
			cacheMsg(nt, best, C_VPASS, "\r");
			cacheMsg(nt, best, (gotBest) ? C_VPASS : C_VITER, "\n");
		} while (delta < 0 || rnd < Boltz_prob);
		//cacheMsg(nt, best, C_VPASS, "\n");

		cacheLock(best);
		if (cacheRead(nt, best) < 0)
			cacheWrite(nt, best);
		cacheUnlock(best);
		cacheRestore(nt, best);
		cacheSave(nt, cacheCost(best), best);

		if (pooling) {
			cacheReset(pool);
			cacheSave(nt, cacheCost(best), pool);
			cacheOpen(pool, "SOLN");
			cacheLock(pool);
			if (cacheRead(nt, pool) < 0)
				cacheWrite(nt, pool);
			cacheUnlock(pool);
			cacheClose(pool);
			cacheSave(nt, cacheCost(pool), best);
		}
		temperature *= 0.95;
	} while (temperature > 0);
	//cacheMsg(nt, best, C_VPASS, "\n");

	cacheClose(best);

	cacheMsg(nt, best, C_VPASS, "%2s ", tag);
	cacheMsg(nt, best, C_VPASS, "T%6d/%6d %C-%R u%.3f ",
		0, 0, cacheCost(best), cacheRootCost(best),
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VPASS, "x%3d     ", cacheMixedNodes(best));
	cacheMsg(nt, best, C_VPASS, "    |: R=%ld, ", RetCost);

	cacheOpen(best, "SOLN");
	cacheLock(best);
	if (cacheRead(nt, best) < 0)
		cacheWrite(nt, best);
	cacheUnlock(best);
	cacheClose(best);
	cacheMsg(nt, best, C_VRUN, "soln: %C-%R ",
		cacheCost(best), cacheRootCost(best));
	cacheMsg(nt, best, C_VRUN, "u%.3f ",
		(double) cacheCost(best)/cacheNodes(best));
	cacheMsg(nt, best, C_VRUN, "x%3d", cacheMixedNodes(best));
	cacheMsg(nt, best, C_VRUN, "\n");

	cacheFree(nt, curr);
	cacheFree(nt, ratchet);
	return cacheCost(best);
}
