
/*
	hyp.c - routines for handling hypothetical characters
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
#include "state.h"
#include "boot.h"

#define FIRST_FOUND 1

static int		hypAncestor(Net *nt, Length basecost, Length rootCost,
					Cursor tt, Cache *best);
static int		hypPoly(Net *nt, Length basecost, Length rootCost,
					Cursor src, Cache *best);
static int		hypLink(Net *nt, Length basecost, Cursor dst, Cache *best);
static int		hypRoot(Net *nt, Cache *best);

static int		hypReplace(Net *nt, Length currcost, Cache *best);

//
// hypAncestor -- resolve direct ancestories to
//	a hypothetical ancestor and cost.
//
//	e.g.  A -> B -> C @ B ==> A -> [1], [1] -> B, [1] -> C
//

int
	hypAncestor(Net *nt, Length basecost, Length rootCost,
		Cursor tt, Cache *best)
{
	Taxa *tx = nt->taxa;
	Length got;
	Cursor hyp;
	Link link[1];
	int cached = NO;
	Cursor t;

	assert( tt < tx->nExtant );

	hyp = hypNew(nt);

	link->from = hyp;
	link->to = tt;
	ntConnect(nt, link);

	// Do reconnections. 
	// O:: Is use of nt->Ups[] and nt->Dns[] faster?
	//     Or could I get into problems if ntDisconnect() changes them?
	for (t = 0; t < nt->maxTax; t++) {
		if (t == hyp)
			continue;
		if (!nt->inuse[t])
			continue;

		if (nt->connection[t][tt]
#if DO_CLINK
		&& t != txCorrecting(tx,tt)
#endif
		) {
			link->from = t;
			link->to = tt;
			ntDisconnect(nt, link);
			link->to = hyp;
			ntConnect(nt, link);
		}

		if (nt->connection[tt][t]
#if DO_CLINK
		&& tt != txCorrecting(tx,t)
#endif
		) {
			link->from = tt;
			link->to = t;
			ntDisconnect(nt, link);
			link->from = hyp;
			ntConnect(nt, link);
		}
	}

	got = basecost;
	if (ntHypConstrained(nt)) {
		got = -1;
		goto break_out;
	}

	if (nt->nParents[tt] == 0)
		got += nt->cumes[tt];
	got = stAncCost(nt, hyp, nt->cumes, got, cacheCost(best));
	if (!cacheBetter(nt, got, rootCost, nt->nLinks, best))
		goto break_out;

	stFixNode(nt, hyp);

#if DO_POLE
	// Pole check all of hyp's kids now
	if (!ntKidsPoleCheck(nt, hyp) || !ntPoleCheck(nt,hyp))
		goto break_out;
#endif

	cached = cacheSave(nt, got, best);
	cacheMsg(nt, best, C_VCACHE, "(%N)=%N %C-%R ", tt, hyp, got, cacheRootCost(best));

break_out:
	// Undo connections.
	// O:: Is use of nt->Ups[] and nt->Dns[] faster?
	for (t = 0; t < nt->maxTax; t++) {
		if (t == tt)
			continue;
		if (!nt->inuse[t])
			continue;
		if (nt->connection[t][hyp]) {
			link->from = t;
			link->to = hyp;
			ntDisconnect(nt, link);
			link->to = tt;
			ntConnect(nt, link);
		}

		if (nt->connection[hyp][t]) {
			link->from = hyp;
			link->to = t;
			ntDisconnect(nt, link);
			link->from = tt;
			ntConnect(nt, link);
		}
	}

	link->from = hyp;
	link->to = tt;
	ntDisconnect(nt, link);
	hypClear(nt, hyp);

	return cached;
}

//
// hypPoly - Resolve an n-polychotomy into an (n-1)-polychotomy
//
//	e.g.  A->B, A->C ==> A->[1], [1]->B, [1]->C
//

int
	hypPoly(Net *nt, Length basecost, Length rootCost, Cursor tt, Cache *best)
{
	Taxa *tx = nt->taxa;
	Length got;
	Cursor hyp;
	Link link[1];
	int cached = NO;
	static Cursor *kids;
	int nKids;
	Cursor k1, k2;
	Cursor t1, t2;
	Cursor dn;
	int nMixed = 0, nUnmixed = 0;
	Cursor root = ntRoot(nt);
	Length newRootCost = 0;
	static vunit *rootBase = 0;
	static vunit *rootRdgs = 0;
	CodeID rootCode = 0;

	NEWMEM(kids, tx->nTotal);
	NEWMEM(rootBase, tx->nVunits);
	NEWMEM(rootRdgs, tx->nVunits);

#if ROOTFIX
	// Save off root info
	if (tt == root) {
		rootCode = nt->codes[root];
		ACPY(rootBase, txBase(tx,root), tx->nVunits);
		ACPY(rootRdgs, txRdgs(tx,root), tx->nVunits);
	}
#endif

	hyp = hypNew(nt);
	assert( hyp != TXNOT );
	nt->cumes[hyp] = 0;

	nKids = 0;
	for (dn = nt->Dns[tt]; dn != -1; dn = nt->br[dn].nxtDn) {
		Cursor to = nt->br[dn].to;
		kids[nKids++] = to;
		
		if (nt->nParents[to] > 1)
			nMixed++;
		else if (!txFrag(tx,to))
			nUnmixed++;
	}
	if (nt->nParents[tt] == 0)
		--nUnmixed;
	if (nt->nParents[tt] > 1)
		nMixed++;

	link->from = tt;
	link->to = hyp;
	ntConnect(nt, link);

	for (k1 = 0; k1 < nKids; k1++) {
		Flag mixed1, frag1;
		t1 = kids[k1];

		mixed1 = (nt->nParents[t1] - 1);
		frag1 = txFrag(tx,t1);

#if DO_CLINK
			if (txCorrecting(tx,t1) == tt)
				continue;
#endif

		link->from = tt;
		link->to = t1;
		ntDisconnect(nt, link);
		link->from = hyp;
		ntConnect(nt, link);

		for (k2 = k1+1; k2 < nKids; k2++) {
			Flag mixed2, frag2;
			int justCached = NO;
			t2 = kids[k2];

			mixed2 = (nt->nParents[t2] - 1);
			frag2 = txFrag(tx,t2);

			// Incremental ntHypConstrained check here:
			if (mixed1 && mixed2) {
				continue;
			} else if (mixed1) {
				if (frag2 || mixed1 > 1) continue;
			} else if (mixed2) {
				if (frag1 || mixed2 > 1) continue;
			} else {
				if (!frag1 && !frag2 && nMixed >= nUnmixed) continue;
			}
#if DO_CLINK
			if (txCorrecting(tx,t2) == tt)
				continue;
#endif

			link->from = tt;
			link->to = t2;
			ntDisconnect(nt, link);
			link->from = hyp;
			ntConnect(nt, link);

			got = stPolyCost(nt, hyp, tt, t1, t2, nt->cumes,
				basecost, cacheCost(best));
#if ROOTFIX
			newRootCost = stPolyRootCost(nt, hyp, tt, root, rootRdgs, rootCost);
#else
			newRootCost = rootCost;
#endif
			if (cacheBetter(nt, got, newRootCost, nt->nLinks, best)) {
				stFixNode(nt, hyp);
				if (newRootCost != rootCost) {
					cacheMsg(nt, best, C_VCACHE, "RC-%d ", rootCost-newRootCost);
					stFixNode(nt, root);
				}
#if DO_POLE
				if (ntKidsPoleCheck(nt, hyp))
#endif
				justCached = cached = cacheSave(nt, got, best);
				if (justCached) {
					cacheMsg(nt, best, C_VCACHE, "%N(%N,%N)=%N %C-%R ",
						tt, t1, t2, hyp, got, cacheRootCost(best));
				}

				// assert( !ntHypConstrained(nt) );

				// Restore old ROOT
				if (newRootCost != rootCost) {
					if (justCached)
						assert( newRootCost == cacheRootCost(best) );
					nt->codes[root] = rootCode;
					ACPY(txBase(tx,root), rootBase, tx->nVunits);
					txPermuteTaxon(tx, root);
				}
			}

			link->to = t2;
			link->from = hyp;
			ntDisconnect(nt, link);
			link->from = tt;
			ntConnect(nt, link);
		}

		link->to = t1;
		link->from = hyp;
		ntDisconnect(nt, link);
		link->from = tt;
		ntConnect(nt, link);
	}
	link->to = hyp;
	ntDisconnect(nt, link);
	hypClear(nt, hyp);

	return cached;
}

int
	hypLink(Net *nt, Length basecost, Cursor dst, Cache *best)
{
	Length got;
	Link link[1];
	int cached = NO;
	Length bestCost = cacheCost(best);

	assert( dst != TXNOT );
	link->from = nt->taxa->nTotal;
	link->to   = nt->taxa->nTotal;

#if DO_MAXMIX
		if (nt->nMixed >= nt->maxMix)
			return cached;
#endif

	got = ntIncLink(nt, basecost, bestCost, dst, link);
	if (got < bestCost) {
		ntConnect(nt, link);
		if (cacheSave(nt, got, best)) {
			cacheMsg(nt, best, C_VCACHE, "%L %C ", link, got);
			cached = YES;
		}
		ntDisconnect(nt, link);
	}
	return cached;
}

int
	hypRoot(Net *nt, Cache *best)
{
	Length got;
	Link link[1];
	int cached = NO;

	link->from = hypNew(nt);
	for (link->to = 0; link->to < link->from; link->to++) {
		if (nt->nParents[link->to] == 0)
			ntConnect(nt, link);
	}
	stFixNode(nt, link->from);
	got = ntCost(nt);
	if (cacheSave(nt, got, best)) {
		cacheMsg(nt, best, C_VSWAP, "ROOT:%N %C ", link->from, got);
		cached = YES;
	}

	for (link->to = 0; link->to < link->from; link->to++) {
		if (nt->connection[link->from][link->to])
			ntDisconnect(nt, link);
	}
	return cached;
}

void
	hypFix(Net *nt)
{
	Cursor *nwrk=0, *node=0, *nend=0;

	nend = ntPostorderNodes(nt, &nwrk);
	for (node = nwrk; node < nend; node++)
		stFixNode(nt, *node);
	free(nwrk);
}

//	hypReplace - try the alternatives, caching the best

static int
	hypReplace(Net *nt, Length basecost, Cache *best)
{
	int cached = NO;
	Cursor tt, root = ntRoot(nt);
	Length rootCost = txApographic(nt->taxa, root, nt->outgroup);

	for (tt = 0; tt < nt->taxa->nExtant; tt++) {
		if (nt->nChildren[tt] < 1)
			continue;
		if (hypAncestor(nt, basecost, rootCost, tt, best))
			cached = YES;
	}

	// Try all polys
	for (tt = 0; tt < nt->maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;				// Taxon not in use (e.g. old hyp)
		if (nt->nChildren[tt] < 2)
			continue;

		if (hypPoly(nt, basecost, rootCost, tt, best))
			cached = YES;
	}

	// Try all links
	for (tt = 0; tt < nt->maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;				// Taxon not in use (e.g. old hyp)
		if (tt == root)
			continue;				// Root is always constrained
		if (nt->nParents[tt] > 0 && nt->cumes[tt] <= RetCost)
			continue;

		if (hypLink(nt, basecost, tt, best))
			cached = YES;
	}

	return cached;
}

static Cursor
	hypGetNew(Net *nt, Cursor hyp)
{
	Cursor h;

	for (h = hyp; h < nt->taxa->nTotal; h++) {
		if (!nt->inuse[h]) {
			nt->inuse[h] = YES;
			if (h >= nt->maxTax)
				nt->maxTax = h+1;
			nt->nTaxa++;
			return h;
		}
	}
	abort();
	return TXNOT;
}

void
	hypClear(Net *nt, Cursor hyp)
{
	if (!nt->inuse[hyp])
		return;
	nt->inuse[hyp] = NO;
	--nt->nTaxa;
	while (!nt->inuse[nt->maxTax-1])
		--nt->maxTax;
}


Cursor
	hypNew(Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor hyp, tt;

	hyp = hypGetNew(nt, tx->nExtant);
	if (hyp == TXNOT)
		return hyp;

	for (tt = tx->nExtant; tt < nt->maxTax; tt++)
		nt->noanc[tt][hyp] = NO;
	ZERO(nt->noanc[hyp],tx->nTotal);
	nt->noanc[hyp][hyp] = YES;

	return hyp;
}

//////////////////////////
// 
// Optimization Routines
//
// CRR - Connection removal & replacement.
//
// Remove a connection (including the hypotheticals nodes) and
//	try a replacement.  Cache the best one.

static Length hypRemAnc(Net *nt, Cache *start, Cache *best);
static Length hypRemPoly(Net *nt, Cache *start, Cache *best);
static Length hypRemLink(Net *nt, Cache *start, Cache *best);

typedef struct RR {
	Length *cumes;
	Flag **vc;
	Flag **vcAlias;
} RR;

static RR *hypRR(Net *nt);
static Length hypAncRR( Net *nt, Cache *start, Cache *best, RR *rr);
static Length hypPolyRR(Net *nt, Cache *start, Cache *best, RR *rr);
static Length hypLinkRR(Net *nt, Cache *start, Cache *best, RR *rr);

Length
	hypCRR(Net *nt, Cache *start, Cache *best, char *prefix)
{
	static Cache *back = 0;
	static Net *ntBack = 0;
	RR *rr;

	Length got, gotBack, got2;

	if (ntBack != nt || !back) {
		back = cacheNew(nt);
		ntBack = nt;
	}

	// If best is new or reset, start it off with start
	cacheRestore(nt, start);
	cacheSave(nt, cacheCost(start), best);

	cacheSetVerbosity(back, cacheVerbosity(best));

	do {
		cacheMark(best);
		//assert( !ntHypConstrained(nt) );
		got = ntCost(nt);	  		// Calc. marginal costs
		assert( got == cacheCost(start) );
		ntPropagate(nt);

		// MORE: Add more links if we can, /*but not if we're perturbed*/
		if (nt->nLinks < 2*(nt->taxa->nExtant - 1)) {
			cacheMsg(nt, best, C_VSWAP, "%s%d/%C-%R MORE:",
				prefix, nt->nLinks, got, cacheRootCost(best));
			hypReplace(nt, got, best);
			cacheMsg(nt, best, C_VSWAP, "\n");
			if (cacheCached(best))
				goto first_found;
		}

		// BACK: Remove the worst link and try again.
		cacheReset(back);
		cacheMsg(nt, best, C_VSWAP, "%s%d/%C-%R BACK:",
			prefix, nt->nLinks, got, cacheRootCost(best));

		gotBack = hypRemLink(nt, start, back);
		cacheMsg(nt, best, C_VCACHE, "%C-%R ", gotBack, cacheRootCost(back));
		gotBack = hypRemPoly(nt, start, back);
		cacheMsg(nt, best, C_VCACHE, "%C-%R ", gotBack, cacheRootCost(back));
		gotBack = hypRemAnc(nt, start, back);
		cacheMsg(nt, best, C_VCACHE, "%C-%R ", gotBack, cacheRootCost(back));

		if (cacheComp(nt, back, best)) {
			cacheRestore(nt, back);
			got = cacheCost(back);
			// Ensure that we will not cache bounded costs.
			cacheReset(best);
			cacheSave(nt, got, best);
			cacheReset(start);
			cacheSave(nt, got, start);
			cacheMsg(nt, best, C_VCACHE, "DROP");
			cacheMsg(nt, best, C_VSWAP, "\n");
			got++;
			goto first_found;
		}

		if (cacheCached(back)) {
			cacheRestore(nt, back);
			got2 = ntCost(nt);	  		// Calc. marginal costs
			assert( got2 == cacheCost(back) );
			assert( gotBack == cacheCost(back) );
		} else
			gotBack = ntCost(nt); 				// Calc. marginal costs
		ntPropagate(nt);
		hypReplace(nt, gotBack, best);
		cacheMsg(nt, best, C_VSWAP, "\n");
		if (cacheCached(best))
			goto first_found;
		cacheRestore(nt, best);
		cacheSave(nt, cacheCost(best), start);

		got2 = ntCost(nt);	  		// Calc. marginal costs
		assert( got2 == cacheCost(best) );
		ntPropagate(nt);
		rr = hypRR(nt);

		// LINK: try to add a link
		cacheMsg(nt, best, C_VSWAP, "%s%d/%C-%R LINK:",
			prefix, nt->nLinks, got, cacheRootCost(best));
		hypLinkRR(nt, start, best, rr);
		cacheMsg(nt, best, C_VSWAP, "\n");
		if (cacheCached(best))
			goto first_found;

		// POLY: try to resolve a polychotomy
		cacheMsg(nt, best, C_VSWAP, "%s%d/%C-%R POLY:",
			prefix, nt->nLinks, got, cacheRootCost(best));
		hypPolyRR(nt, start, best, rr);
		cacheMsg(nt, best, C_VSWAP, "\n");
		if (cacheCached(best))
			goto first_found;

		// ANCS: try to posit an ancestor
		cacheMsg(nt, best, C_VSWAP, "%s%d/%C-%R ANCS:",
			prefix, nt->nLinks, got, cacheRootCost(best));
		hypAncRR(nt, start, best, rr);
		cacheMsg(nt, best, C_VSWAP, "\n");

	first_found:
		cacheRestore(nt, best);
		cacheSave(nt, cacheCost(best), start);
	} while (cacheCached(best));

	return cacheCost(best);
}

/*
	hypTransfer{Up,Down} - code common to hyp{Anc,Poly}RR and hypRem{Anc,Poly}
		for transfering parents and kids of kid from node 'tt'
		to 't1'.  Returns whether to "break out."

	hypTransferUp   - assumes t1->t2, so ntConstrained(nt, tt, to) is NO
	hypTransferDown - assumes tt->t1, so ntConstrained(nt, fr, t1) is NO

	hypUntransfer - reverses hypTransfer.  No return.
*/

static int
	hypTransferUp(Net *nt, Cursor tt, Cursor t1, Link *work, Link *end)
{
	Link *link2;

	for (link2 = work; link2 < end; link2++) {
		// Transfer tt's parents to be t1's parents
		if (link2->to == tt && link2->from != t1) {
			if (ntConstrained(nt, link2->from, t1))
				return YES;
			ntDisconnect(nt, link2);
			link2->to = t1;
			ntConnect(nt, link2);
			link2->to = tt;
		}

		// Transfer tt's children to be t1's children
		if (link2->from == tt && link2->to != t1) {
			// Always not constrained
			ntDisconnect(nt, link2);
			link2->from = t1;
			ntConnect(nt, link2);
			link2->from = tt;
		}
	}
	return NO;
}

static int
	hypTransferDown(Net *nt, Cursor tt, Cursor t1, Link *work, Link *end)
{
	Link *link2;

	for (link2 = work; link2 < end; link2++) {
		// Transfer tt's parents to be t1's parents
		if (link2->to == tt && link2->from != t1) {
			// Always not constrained
			ntDisconnect(nt, link2);
			link2->to = t1;
			ntConnect(nt, link2);
			link2->to = tt;
		}

		// Transfer tt's children to be t1's children
		if (link2->from == tt && link2->to != t1) {
			if (ntConstrained(nt, t1, link2->to))
				return YES;
			ntDisconnect(nt, link2);
			link2->from = t1;
			ntConnect(nt, link2);
			link2->from = tt;
		}
	}
	return NO;
}

static void
	hypUntransfer(Net *nt, Cursor tt, Cursor t1, Link *work, Link *end)
{
	Link *link2;

	// Restore links
	for (link2 = work; link2 < end; link2++) {
		// Untransfer tt's parents to be t1's parents
		if (link2->to == tt) {
			link2->to = t1;
			ntDisconnect(nt, link2);
			link2->to = tt;
			ntConnect(nt, link2);
		}

		// Untransfer tt's children to be t1's children
		if (link2->from == tt && link2->to != t1) {
			link2->from = t1;
			ntDisconnect(nt, link2);
			link2->from = tt;
			ntConnect(nt, link2);
		}
	}
	
	// In case tt and t1 have same parents or kids
	// Mop up any link with t1 (also reconnects link1)
	for (link2 = work; link2 < end; link2++) {
		if (link2->from == t1 || link2->to == t1)
			ntConnect(nt, link2);
	}

}

/*
	hypCostNode, HYP_IncCost, HYP_RestoreIncCost
		... for Incremental Node Cost calculation
*/

static RR *
	hypRR(Net *nt)
{
	Taxa *tx = nt->taxa;
	static RR rr[1];
	static int init = NO;

	if (!init)  {
		const int alignTotal  = (tx->nTotal +ALIGN) & ~ALIGN;
		const int alignVunits = (tx->nVunits+ALIGN) &~ ALIGN;

		newmem(rr->cumes, tx->nTotal);
		NEWMAT(rr->vc, alignTotal, alignVunits);
		ASET(rr->vc[0], 0, alignTotal * alignVunits);
		newmem(rr->vcAlias, tx->nTotal);
		init = YES;
	}
	ACPY(rr->cumes, nt->cumes, tx->nTotal);
	return rr;
}

static Length
	hypCostNode(Net *nt, Taxa *tx, Cursor t, RR *rr)
{
	Cursor up = nt->Ups[t];
	Length got = tx->nVunits;

	rr->vcAlias[t] = nt->vcs[t]; nt->vcs[t] = rr->vc[t];
	if (up == -1) return 0;
	got -= stLinkCost(tx, nt->br[up].fr, t, nt->vcs[t]);
	while ((up = nt->br[up].nxtUp) != -1) {
		got += RetCost;
		got -= stRetLinkCost(tx, nt->br[up].fr, t, nt->vcs[t]);
	}
	return got;
}
#define HYP_IncCost(got, nt, tx, t, rr) do { \
	Cursor dn; \
	got += (nt->cumes[t] = hypCostNode(nt, tx, t, rr)) - rr->cumes[t]; \
	for (dn = nt->Dns[t]; dn != -1; dn = nt->br[dn].nxtDn) { \
		Cursor to = nt->br[dn].to; \
		if (nt->cumes[to] != rr->cumes[to]) continue; \
		got -= rr->cumes[to]; \
		got += (nt->cumes[to] = hypCostNode(nt, tx, to, rr)); \
	} \
} while (0)

#define HYP_RestoreIncCost(nt, tx, t, rr) do { \
	Cursor dn; \
	nt->cumes[t] = rr->cumes[t]; \
	nt->vcs[t] = rr->vcAlias[t]; \
	for (dn = nt->Dns[t]; dn != -1; dn = nt->br[dn].nxtDn) { \
		Cursor to = nt->br[dn].to; \
		nt->cumes[to] = rr->cumes[to]; \
		nt->vcs[to] = rr->vcAlias[to]; \
	} \
} while (0)


static Length
	hypAncRR(Net *nt, Cache *start, Cache *best, RR *rr)
{
	Taxa *tx = nt->taxa;
	int cached = NO;
	Length bestCost = cacheCost(best);
	Link *work, *link, *end;
	Cursor t;
	Length oldBase, newBase;

	oldBase = cacheCost(start);

	// Bring up each descendant
	end = cacheListLinks(start, &work);
	for (link = work; link < end; link++) {
		Cursor tt = link->from;
		Cursor t1 = link->to;
		//int wasRoot;

		if (tt < tx->nExtant || t1 >= tx->nExtant)
			continue;
		
		//wasRoot = (nt->nParents[tt] == 0);

		// Got a descendant: it is 't1'.  Now fix links.
		ntDisconnect(nt, link);
		if (hypTransferDown(nt, tt, t1, work, end))
			goto break_out;

#if DO_POLE
		// Pole check all of t1's kids now
		if (!ntKidsPoleCheck(nt, t1))
			goto break_out;
#endif

		hypClear(nt, tt);
		if (!ntHypConstrained(nt)) {
			newBase = oldBase;
			newBase -= rr->cumes[tt];

			HYP_IncCost(newBase, nt, tx, t1, rr);
			cached = hypAncestor(nt, newBase, -1, t1, best);

			// Be verbose
			if (cached)
				cacheMsg(nt, best, C_VLINK, "{%N<%N} ", tt, t1);

			HYP_RestoreIncCost(nt, tx, t1, rr);
			cacheNodeRestore(nt, t1, start);
		}

		t = hypGetNew(nt, tt);
		assert( t == tt );
		cacheNodeRestore(nt, tt, start);

	break_out:
		hypUntransfer(nt, tt, t1, work, end);

#if FIRST_FOUND
		// Return first found
		if (cached && cacheCost(best) < bestCost)
			break;
#endif
	}

	return cacheCost(best);
}

static Length
	hypRemAnc(Net *nt, Cache *start, Cache *best)
{
	Taxa *tx = nt->taxa;
	int cached = NO;
	Length basecost = cacheCost(start);
	Length rootCost = txApographic(tx, ntRoot(nt), nt->outgroup);
	Link *work, *link, *end;
	Cursor t;

	// Bring up each descendant
	end = cacheListLinks(start, &work);
	for (link = work; link < end; link++) {
		Cursor tt = link->from;
		Cursor t1 = link->to;
		int justCached = NO;
		int wasRoot = (nt->nParents[tt] == 0);

		if (tt < tx->nExtant)
			continue;

		// Got a descendant: it is 't1'.  Now fix links.
		ntDisconnect(nt, link);
		if (hypTransferDown(nt, tt, t1, work, end))
			goto break_out;

		hypClear(nt, tt);
		if (!ntHypConstrained(nt)) {
			Length got = basecost;

			got -= nt->cumes[t1];
			if (!wasRoot)
				got -= nt->cumes[tt];
			got = stAncCost(nt, t1, nt->cumes, got, -1);
			if (cacheBetter(nt, got, (wasRoot) ? -1 : rootCost, nt->nLinks, best)) {
				stFixNode(nt, t1);
#if DO_POLE
				if (ntKidsPoleCheck(nt, t1) && ntPoleCheck(nt, t1))
#endif
				justCached = cached = cacheSave(nt, got, best);
				cacheNodeRestore(nt, t1, start);
			}

			// Be verbose
			if (justCached)
				cacheMsg(nt, best, C_VLINK, "{%N<%N} ", tt, t1);
		}

		t = hypGetNew(nt, tt);
		assert( t == tt );

	break_out:
		hypUntransfer(nt, tt, t1, work, end);
	}

	return cacheCost(best);
}

static Length
	hypPolyRR(Net *nt, Cache *start, Cache *best, RR *rr)
{
	Taxa *tx = nt->taxa;
	Cursor tt, t1;
	int cached = NO;
	Link link[1];
	Length bestCost = cacheCost(best);
	Link *work, *end;
	Cursor t;
	Length oldBase, newBase;

	oldBase = cacheCost(start);

	end = cacheListLinks(start, &work);
	for (tt = tx->nExtant; tt < nt->maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;
		
		// Push up to each parent
		for (t1 = 0; t1 < nt->maxTax; t1++) {
			int justCached = NO;

			if (!nt->connection[t1][tt])
				continue;
			if (!nt->inuse[t1])
				continue;
			
			// Got a parent: it is 't1'.  Now fix links.
			link->from = t1;
			link->to = tt;
			ntDisconnect(nt, link);

			if (hypTransferUp(nt, tt, t1, work, end))
				goto break_out;

			hypClear(nt,tt);
			if (!ntHypConstrained(nt)) {
				stFixNode(nt, t1);
				newBase = oldBase;
				HYP_IncCost(newBase, nt, tx, tt, rr);
				HYP_IncCost(newBase, nt, tx, t1, rr);
#if DO_POLE
				if (ntPoleCheck(nt, t1) && ntKidsPoleCheck(nt, t1))
#endif
				justCached = cached = hypPoly(nt, newBase, -1, t1, best);

				// Be verbose
				if (justCached)
					cacheMsg(nt, best, C_VLINK, "{%N>%N} ", tt, t1);

				HYP_RestoreIncCost(nt, tx, t1, rr);
				HYP_RestoreIncCost(nt, tx, tt, rr);
				cacheNodeRestore(nt, t1, start);
			}

			t = hypGetNew(nt, tt);
			assert( t == tt );
			cacheNodeRestore(nt, tt, start);

		break_out:
			hypUntransfer(nt, tt, t1, work, end);

#if FIRST_FOUND
			// Return first found
			if (cached && cacheCost(best) < bestCost)
				goto first_found;
#endif
		}
	}

first_found:
	return cacheCost(best);
}

static Length
	hypRemPoly(Net *nt, Cache *start, Cache *best)
{
	Taxa *tx = nt->taxa;
	Cursor tt, t1;
	int cached = NO;
	Link link[1];
	Length basecost = cacheCost(start);
	Length rootCost = txApographic(tx, ntRoot(nt), nt->outgroup);
	Link *work, *end;
	Cursor t;

	end = cacheListLinks(start, &work);
	for (tt = tx->nExtant; tt < nt->maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;
		
		// Push up to each parent
		for (t1 = 0; t1 < nt->maxTax; t1++) {
			int justCached = NO;
			int wasRoot = (nt->nParents[t1] == 0);

			if (!nt->connection[t1][tt])
				continue;
			if (!nt->inuse[t1])
				continue;
			
			// Got a parent: it is 't1'.  Now fix links.
			link->from = t1;
			link->to = tt;
			ntDisconnect(nt, link);

			if (hypTransferUp(nt, tt, t1, work, end))
				goto break_out;

			hypClear(nt,tt);
			if (!ntHypConstrained(nt)) {
				Length got = basecost;
				if (!wasRoot)
					got -= nt->cumes[t1];
				got -= nt->cumes[tt];
				got = stAncCost(nt, t1, nt->cumes, got, cacheCost(best));
				if (cacheBetter(nt, got, (wasRoot) ? -1 : rootCost, nt->nLinks, best)) {
					stFixNode(nt, t1);
#if DO_POLE
					if (ntKidsPoleCheck(nt, t1) && ntPoleCheck(nt, t1))
#endif
					justCached = cached = cacheSave(nt, got, best);
					cacheNodeRestore(nt, t1, start);
				}

				// Be verbose
				if (justCached)
					cacheMsg(nt, best, C_VLINK, "{%N>%N} ", tt, t1);
			}

			t = hypGetNew(nt, tt);
			assert( t == tt );

		break_out:
			hypUntransfer(nt, tt, t1, work, end);
		}
	}

	return cacheCost(best);
}

static Length
	hypLinkRR(Net *nt, Cache *start, Cache *best, RR *rr)
{
	Taxa *tx = nt->taxa;
	int cached = NO;
	Length bestCost = cacheCost(best);
	Length oldBase, newBase;
	Branch *br, *brEnd = &nt->br[2*tx->nExtant];

	oldBase = cacheCost(start);

	for (br = nt->br; br < brEnd; br++) {
		Cursor t1, t2;
		Link link[1];
		int justCached;

		if (br->nxtBr != -2)
			continue;

		t1 = br->fr;
		t2 = br->to;
		
		if (nt->nParents[t2] == 1)
			continue;

#if DO_CLINK
		if (txCorrecting(tx,t2) == t1)
			continue;
#endif

		justCached = NO;
		link->from = t1;
		link->to = t2;
		ntDisconnect(nt, link);

		stFixNode(nt, t2);
		stFixNode(nt, t1);

		newBase = oldBase;

		HYP_IncCost(newBase, nt, tx, t1, rr);
		HYP_IncCost(newBase, nt, tx, t2, rr);

#if DO_POLE
		// Pole check all of t1's kids now
		if (ntKidsPoleCheck(nt, t1))
#endif

		if (nt->cumes[t2] > RetCost && hypLink(nt, newBase, t2, best))
			justCached = cached = YES;

		// Be verbose
		if (justCached)
			cacheMsg(nt, best, C_VLINK, "{%N:>%N} ", t1, t2);

		HYP_RestoreIncCost(nt, tx, t1, rr);
		HYP_RestoreIncCost(nt, tx, t2, rr);

		cacheNodeRestore(nt, t1, start);
		cacheNodeRestore(nt, t2, start);
		ntConnect(nt, link);

#if FIRST_FOUND
		// Return first found
		if (cached && cacheCost(best) < bestCost)
			goto first_found;
#endif
	}

first_found:
	return cacheCost(best);
}

static Length
	hypRemLink(Net *nt, Cache *start, Cache *best)
{
	Taxa *tx = nt->taxa;
	Cursor t1, t2;
	int cached = NO;
	Length basecost = cacheCost(start);
	Length rootCost = txApographic(tx, ntRoot(nt), nt->outgroup);
	static Length *cumes;

	if (!cumes) {
		newmem(cumes, tx->nTotal);
	}

	ACPY(cumes, nt->cumes, tx->nTotal);
	for (t1 = 0; t1 < nt->maxTax; t1++) {
		int wasRoot = (nt->nParents[t1] == 0);
		if (!nt->inuse[t1])
			continue;
		
		for (t2 = 0; t2 < nt->maxTax; t2++) {
			int justCached = NO;
			Link link[1];
			Length got = basecost;

			if (!nt->connection[t1][t2])
				continue;
			if (!nt->inuse[t2])
				continue;
			if (nt->nParents[t2] == 1)
				continue;
			
#if DO_CLINK
			if (txCorrecting(tx,t2) == t1)
				continue;
#endif

			link->from = t1;
			link->to = t2;
			ntDisconnect(nt, link);

			got -= cumes[t2];
			if (nt->nParents[t2] == 0)
				got += tx->nVunits;
			got = stAncCost(nt, t2, cumes, got, -1);
			if (!wasRoot)
				got -= nt->cumes[t1];
			got = stAncCost(nt, t1, cumes, got, -1);

			if (cacheBetter(nt, got, (wasRoot) ? -1 : rootCost, nt->nLinks, best)) {
				stFixNode(nt, t1);
				stFixNode(nt, t2);
#if DO_POLE
				if (ntKidsPoleCheck(nt, t2))
#endif
				justCached = cached = cacheSave(nt, got, best);
				cacheNodeRestore(nt, t1, start);
				cacheNodeRestore(nt, t2, start);
			}

			// Be verbose
			if (justCached) {
				cacheMsg(nt, best, C_VLINK, "{%N:>%N} ", t1, t2);
			}

			ntConnect(nt, link);
		}
	}

	return cacheCost(best);
}

extern Length
	hypLinkUp(Net *nt, Cache *curr, Cache *best)
{
	int n;
	int nLeft = nt->taxa->nExtant - 1;	// Number of taxa left to link.

#if DO_CLINK
	// Link up correctors to their correcting manuscript
	for (n = 0; n < nt->taxa->nExtant; n++) {
		Link link[1];
		link->from = txCorrecting(nt->taxa,n);
		link->to = n;

		if (link->from == TXNOT)
			continue;
		ntConnect(nt, link);

		Length got = ntCost(nt);
		cacheSave(nt, got, curr);
		cacheMsg(nt, curr, C_VSWAP, "%N->%N/%C\n", link->from, link->to, got);
	}
#endif

	// Initial link up
	for (n = 0; n < nLeft; n++) {
		Length got, cost;
		Cursor to;

		cost = ntCost(nt);		// Needed to set up incrementals
		ntPropagate(nt);
		assert( cost == cacheCost(curr) );

		cacheMsg(nt, best, C_VSWAP, "%d/%C ", nt->nLinks, cost);
		for (to = 0; to < nt->maxTax; to++) {
			if (nt->nLinks < nt->taxa->nExtant-1 && nt->nParents[to] > 0)
				continue;				// Initial link up.
			if (!nt->inuse[to])
				continue;				// Taxon not in use (e.g. old hyp)

			hypLink(nt, cost, to, best);
		}

		got = cacheRestore(nt, best);

		// Couldn't add, so just link all the orphans up.
		if (got >= cost) {
			hypRoot(nt, best);
			got = cacheRestore(nt, best);
			cacheSave(nt, got, curr);
			cacheMsg(nt, best, C_VSWAP, "\n");
			break;
		}
		cacheMsg(nt, best, C_VSWAP, "\n");

		cost = got;
		cacheSave(nt, got, curr);
	}
	return cacheCost(best);
}
