
/*
	boot.c - routines for statistical analysis (boot strap)
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "netmsg.h"
#include "cache.h"
#include "hyp.h"
#include "boot.h"

typedef struct pt_trie PtTrie;
static PtTrie *	ptNew();
static long		ptTally(PtTrie *pt, Flag *tset, size_t len, int inc);
static void		ptWalk(PtTrie *pt, Flag *tset, Flag *curset, size_t len,
					void *ctx, void (*fcn)(void *, long count));
static void		ptFree(PtTrie *pt);

typedef struct {
	Cursor node;		// Node in the solution
	Flag inuse;			// In use in the solution
	Flag dupTset;		// Duplicates another's tset
	Flag *tset;			// Terminal set
	long mixed;			// No. of times node is mixed
	long *kids;			// Count of direct descendants
} BtNode;

typedef struct {
	FILE *log;			// File to print to
	Flag *curset;		// Current set
	long thresh;		// Threshold
	long reps;			// Replicates
	Cursor len;			// Length of set
	Net *nt;			// Network solution
	Cursor node;		// Current Node
} BtPrintSub;

struct Boot {
	Cache *soln;		// Cached solution
	long reps;			// Number of replicates in boot strap
	long tsets;			// Number of tsets found
	BtNode *nodes;		// Solution nodes
	PtTrie *pt;			// Partition Trie to hold all tsets
};


/*
	bootSubset(a, b, n) - return YES if a <= b, else NO
*/
static int
	bootSubset(Flag *a, Flag *b, int n)
{
	assert( n > 0 );
	while (*a++ <= *b++ && --n > 0)
		;
	return (n == 0);
}

/*
	bootDupTset(nt, node) - return YES if node's tset duplicates another.
		Determined by seeing if the mixed kids is included in the other.
*/
static int
	bootDupTset(Net *nt, Cursor node)
{
	Taxa *tx = nt->taxa;
	Cursor first, second;			// First and second kids
	Cursor dn;

	// Must have two kids
	if (nt->nChildren[node] != 2)
		return NO;
	
	dn = nt->Dns[node];
	first = nt->br[dn].to;
	dn = nt->br[dn].nxtDn;
	second = nt->br[dn].to;

	if (nt->nParents[first] > 1
	&& bootSubset(nt->descendents[first], nt->descendents[second], tx->nExtant))
		return YES;
	if (nt->nParents[second] > 1
	&& bootSubset(nt->descendents[second], nt->descendents[first], tx->nExtant))
		return YES;
	return NO;
}

/*
	bt = bootNew(nt, soln) - Allocate bt for a cached solution (soln).
*/
Boot *
	bootNew(Net *nt, Cache *soln)
{
	Taxa *tx;
	Cursor tt;
	Boot *bt;

	cacheRestore(nt, soln);
	ntPropagate(nt);
	tx = nt->taxa;

	newmem(bt, 1);
	bt->soln = soln;
	bt->reps = 0;
	bt->tsets = 0;

	newmem(bt->nodes, tx->nTotal);
	for (tt = 0; tt < tx->nTotal; tt++) {
		BtNode *bn = &bt->nodes[tt];
		bn->node = tt;
		bn->inuse = nt->inuse[tt];
		bn->dupTset = bootDupTset(nt, tt);
		newmem(bn->tset, tx->nExtant);
		memcpy(bn->tset, nt->descendents[tt], tx->nExtant);
		bn->mixed = 0;
		newmem(bn->kids, tx->nExtant);
		ZERO(bn->kids, tx->nExtant);
	}

	bt->pt = ptNew();
	return bt;
}

/*
	bootCount(nt, bt) - Count BootStrap stats:
		count = no. of times that the node's tset is found in replicants.
		mixed = no. of times that a leaf is mixed.
		kids  = no. of times that a leaf has kids (for each kid)
*/
void
	bootCount(Net *nt, Boot *bt)
{
	Taxa *tx = nt->taxa;
	Cursor tt;

	bt->reps++;
	ntPropagate(nt);

	// Go through each node in the replicate tree
	for (tt = 0; tt < nt->maxTax; tt++) {
		BtNode *leaf;

		if (!nt->inuse[tt])
			continue;
		if (bootDupTset(nt, tt))
			continue;

		if (ptTally(bt->pt, nt->descendents[tt], tx->nExtant, 1) == 1)
			bt->tsets++;

		if (tt >= tx->nExtant)
			continue;
		leaf = &bt->nodes[tt];

		// Count if Mixed
		if (nt->nParents[tt] > 1)
			leaf->mixed++;

		// Count kids of Solution Leaf Nodes
		if (nt->nChildren[tt] > 0) {
			Cursor kk;
			Flag *dd = nt->descendents[tt];
			for (kk = 0; kk < tx->nExtant; kk++) {
				assert( nt->inuse[kk] );
				if (kk == tt)
					continue;
				leaf->kids[kk] += dd[kk];
			}
		}
	}
}

/*
	bootStrap(nt, bt, nReps) - Perform the BootStrap for nReps replicants
*/
void
	bootStrap(Net *nt, Boot *bt, int nReps)
{
	static char prefix[80];
	int nR;							// Rep index
	Cache *curr;					// Resampled start
	Cache *rep;						// Replicate
	Length got;	

	// Set up cache for replicate solution
	rep = cacheNew(nt);
	curr = cacheNew(nt);
	cacheSetVerbosity(rep, C_VITER);

	// Run the Boot strap
	for (nR = 0; nR < nReps; nR++) {
		cacheMsg(nt, bt->soln, C_VPASS, "\rBootstrap replication: %d of %d",
			nR+1, nReps);

		cacheRestore(nt, bt->soln);
		cacheReset(rep);
		cacheReset(curr);
		txPerturb(nt->taxa, nt->codes);
		hypFix(nt);
		got = ntCost(nt);
		cacheSave(nt, got, curr);
		cacheSave(nt, got, rep);
		sprintf(prefix, "Boot[%04d/%04d] ", nR+1, nReps);
		hypCRR(nt, curr, rep, prefix);
		cacheMsg(nt, rep, C_VSWAP, "\n");
		
		// Collect stats
		cacheRestore(nt, rep);
		bootCount(nt, bt);
	}
	cacheMsg(nt, bt->soln, C_VPASS, "\n");

	cacheFree(nt, curr);
	cacheFree(nt, rep);
	return;
}

static void
	bootPrintTset(FILE *log, Taxa *tx, Flag *tset)
{
	Cursor tt;

	for (tt = 0; tt < tx->nExtant; tt++) {
		if (tset[tt])
			fprintf(log, " "NAM_F"", txName(tx,tt));
	}
}

static void
	bootPrintSub(void *ctx, long count)
{
	BtPrintSub *bp = ctx;
	Net *nt = bp->nt;
	Taxa *tx = nt->taxa;
	Cursor node = bp->node, dn;

	if (count < bp->thresh)
		return;

	// Only print subsets that cannot be contained in a child
	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn) {
		Cursor kk = nt->br[dn].to;
		if (bootSubset(bp->curset, nt->descendents[kk], bp->len))
			return;
	}

	// Don't print if curset == tset (i.e. if tset <= curset)
	if (bootSubset(nt->descendents[node], bp->curset, bp->len))
		return;

	fprintf(bp->log, "    subset: % 6ld %4.2f",
		count, (double) count / bp->reps);
	bootPrintTset(bp->log, tx, bp->curset);
	fprintf(bp->log, "\n");
}

/*
	bootPrint(log, nt, bt) - Print boot stats to log
*/
void	
	bootPrint(FILE *log, Net *nt, Boot *bt)
{
	Cursor tt, t2;
	Taxa *tx = nt->taxa;
	BtPrintSub bp;

	bp.log = log;
	newmem(bp.curset, tx->nExtant);
	bp.thresh = 0.007 * bt->reps;
	bp.reps = bt->reps;
	bp.len = tx->nExtant;
	bp.nt = nt;

	cacheRestore(nt, bt->soln);
	ntPropagate(nt);

	fprintf(log, "BOOTSTRAP RESULTS:\n");
	fprintf(log, "reps=%ld\n", bt->reps);
	fprintf(log, "tsets=%ld\n", bt->tsets);

	for (tt = 0; tt < nt->maxTax; tt++) {
		BtNode *bn;
		if (!nt->inuse[tt])
			continue;

		bn = &bt->nodes[tt];
		fprintf(log, "\nNode: "NAM_F"\n", txName(tx,bn->node));

		// Print Terminal Sets:
		if (bn->dupTset) {
			fprintf(log, "      tset:");
			bootPrintTset(log, tx, bn->tset);
			fprintf(log, "\n");
			fprintf(log, "Duplicate terminal set.\n");
		} else {
			long count = ptTally(bt->pt, bn->tset, tx->nExtant, 0);

			fprintf(log, "     count: % 6ld %4.2f",
				count, (double) count / bt->reps);
			bootPrintTset(log, tx, bn->tset);
			fprintf(log, "\n");

			// Print subsets
			bp.node = bn->node;
			ptWalk(bt->pt, bn->tset, bp.curset, tx->nExtant, &bp, bootPrintSub);
		}

		if (tt >= tx->nExtant)
			continue;

		// Print leaf node stats
		fprintf(log, "      leaf: mixed=% 6ld (%4.2f); ",
			bn->mixed, (double) bn->mixed / bt->reps);

		// Print Kid Stats:
		fprintf(log, " kids:");
		for (t2 = 0; t2 < tx->nExtant; t2++) {
			if (bn->kids[t2] <= bp.thresh)
				continue;
			fprintf(log, " "NAM_F"=%ld", txName(tx,t2), bn->kids[t2]);
		}
		fprintf(log, "\n");
	}

	free(bp.curset);
	return;
}

void
	bootFree(Net *nt, Boot *bt)
{
	Cursor tt;
	for (tt = 0; tt < nt->taxa->nTotal; tt++)
		free(bt->nodes[tt].tset);
	free(bt->nodes);
	ptFree(bt->pt);
	free(bt);
}

/////////////////////////////////////////////////

/* Record all partitions */

/*
	Data structure: trie so the tset is not directly
	encoded.  Why?  O(lg N) look up times.  Tset are not
	compact.  Compacting and/or hash requires pass over
	entire string any way.

	Issue: best chunk size?  Best is 1.
*/

#define UNSIZ 1024

union unode {
	union unode *link;
	long leaf;
};
typedef union unode Unode;

struct pool {
	struct pool *link;			// Chain these up
	Unode *next;				// Next one to allocate
	Unode unodes[UNSIZ];		// Pool of unodes
};
typedef struct pool Pool;

struct pt_trie {
	Unode *root;
	Pool *pool;
};

/*
	pl = poolNew(pt) - Allocate and push new pool.
*/
static Pool *
	poolNew(Pool *link)
{
	Pool *pl;

	newmem(pl, 1);
	pl->link = link;
	pl->next = pl->unodes;
	return pl;
}

/*
	poolFree(pt) - Pop and free pool.
*/
static Pool *
	poolFree(Pool *pl)
{
	Pool *link;
	if (!pl)
		return pl;
	link = pl->link;
	free(pl);
	return link;
}
	
/*
	unodeNew(pt, n, is_link) - allocate n unodes from pools, and
		initialize based on whether is link or leaf.
*/
static Unode *
	unodeNew(PtTrie *pt, size_t n, int is_link)
{
	Pool *pl = pt->pool;
	Unode *un;

	if (pl->next - pl->unodes + n > UNSIZ) {
		// Pool size should have been an exact multiple.
		assert( pl->next - pl->unodes == UNSIZ );
		pl = pt->pool = poolNew(pl);
	}
	un = pl->next;
	pl->next += n;

	// Too bad I can't assume that (Unode *) 0 == (long) 0!
	if (is_link) {
		do {
			un[--n].link = 0;
		} while (n > 0);
	} else {
		do {
			un[--n].leaf = 0;
		} while (n > 0);
	}
	return un;
}

/*
	ptTally(pt, tset, len, inc) - increment tset by inc and
		return result.  If inc==0, don't allocate storage,
		so it can be used to fetch values out of sparse set.
*/
static long
	ptTally(PtTrie *pt, Flag *tset, size_t len, int inc)
{
	Unode *un = pt->root;
	size_t off = 0;
	
	while (off = *tset++, --len > 0) {
		// Allocate Node only when needed.
		if (un[off].link == 0) {
			if (inc == 0)
				return 0;
			un[off].link = unodeNew(pt, 1 << 1, len > 0);
		}
		un = un[off].link;
	}

	// At leaf node
	return un[off].leaf += inc;
}

/*
	ptNew() - Allocate new Partition Trie.
*/
static PtTrie *
	ptNew()
{
	PtTrie *pt;

	newmem(pt, 1);
	pt->pool = 0;
	pt->pool = poolNew(pt->pool);
	pt->root = unodeNew(pt, 1 << 1, YES);

	return pt;
}

/*
	ptWalk0 - recursive walking function, assumes len >= 1.
*/
static void	
	ptWalk0(Unode *un, Flag *tset, Flag off, Flag *curset, size_t len,
		void *ctx, void (*fcn)(void *, long count))
{
	// Skip if not a subset
	if (off > *tset)
		return;

	// Copy current element flag to current set
	*curset = off;

	// Base case: at leaf
	if (len == 1) {
		if (un[off].leaf > 0)
			(*fcn)(ctx, un[off].leaf);
		return;
	}

	// Recurse.
	if (un[off].link != 0) {
		ptWalk0(un[off].link, tset+1, 0, curset+1, len-1, ctx, fcn);
		ptWalk0(un[off].link, tset+1, 1, curset+1, len-1, ctx, fcn);
	}
	
}

/*
	ptWalk(pt, tset, len, ctx, fcn) - Walk Partition Trie; visit
		every partition that is a subset of tset.
*/
static void
	ptWalk(PtTrie *pt, Flag *tset, Flag *curset, size_t len, void *ctx,
		void (*fcn)(void *, long count))
{
	ptWalk0(pt->root, tset, 0, curset, len, ctx, fcn);
	ptWalk0(pt->root, tset, 1, curset, len, ctx, fcn);
}

/*
	ptFree(pt) - Free pt
*/
static void
	ptFree(PtTrie *pt)
{
	if (!pt)
		return;
	while (pt->pool)
		pt->pool = poolFree(pt->pool);
	free(pt);
}
