
/*
	net.c - routines for CTA parsimony searching
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
#include "state.h"

#define SORTLINKS 1

//////////////////////////////
//
// ntNew
//

Net *
	ntNew(Taxa *tx)
{
	Net *nt;
	Cursor from;

	// Compute aligned for fostering 32-bit memory moves
	int alignTotal  = (tx->nTotal +ALIGN) & ~ALIGN;
	int alignVunits = (tx->nVunits+ALIGN) &~ ALIGN;

	newmem(nt, 1);

	nt->taxa = tx;
	nt->nLinks = 0;
	nt->nMixed = 0;

	// Some dynamic values related to hypothetical nodes;
	nt->nTaxa = tx->nExtant;
	nt->maxTax = tx->nExtant;
	newmem(nt->inuse, tx->nTotal);
	ASET(nt->inuse, 1, tx->nExtant);
	ASET(nt->inuse+tx->nExtant, 0, tx->nTotal - tx->nExtant);
	nt->outgroup = txFind(tx, getenv("OUTGRP"));

	// Allocate cost and state cache
	NEWMAT(nt->vcBase, alignTotal, alignVunits);
	ASET(nt->vcBase[0], 0, alignTotal * alignVunits);
	newmem(nt->vcs, alignTotal);
	newmem(nt->codes, alignTotal);
	for (from = 0; from < alignTotal; from++) {
		nt->codes[from] = TaxonCode++;
		nt->vcs[from] = nt->vcBase[from];
	}
	newmem(nt->cumes, tx->nTotal);
	ASET(nt->cumes, 0, tx->nTotal);
#if DO_POLE
	newmem(nt->poles, tx->nTotal);
	ASET(nt->poles, 0, tx->nTotal);
	for (Cursor to = 0; to < tx->nExtant; to++)
		nt->poles[to] = txPole(tx, to);
	newmem(nt->banMixed, tx->nTotal);
	ASET(nt->banMixed, 0, tx->nTotal);
#endif

	// Align and linearize connection matrix for ntPropagate()
	NEWMAT(nt->connection, alignTotal, alignTotal);
	ASET(nt->connection[0], 0, alignTotal * alignTotal);

	newmem(nt->br, 2*tx->nExtant);
	newmem(nt->Ups, tx->nTotal);
	newmem(nt->Dns, tx->nTotal);
	for (from = 0; from < 2*tx->nExtant-1; from++)
		nt->br[from].nxtBr = from+1;
	nt->br[from].nxtBr = -1;
	nt->freeBr = 0;
	ASET(nt->Ups, -1, tx->nTotal);
	ASET(nt->Dns, -1, tx->nTotal);

	newmem(nt->nParents, tx->nTotal);
	ASET(nt->nParents, 0, tx->nTotal);
	newmem(nt->nChildren, tx->nTotal);
	ASET(nt->nChildren, 0, tx->nTotal);

	NEWMAT(nt->noanc, alignTotal, alignTotal);
	ASET(nt->noanc[0], 0, alignTotal * alignTotal);

	NEWMAT(nt->descendents, alignTotal, alignTotal);
	ASET(nt->descendents[0], 0, alignTotal * alignTotal);

	for (from = 0; from < tx->nTotal; from++) {
		nt->codes[from] = TaxonCode++;				// R:: redundant?
		nt->noanc[from][from] = YES;

		nt->vcs[from] = nt->vcBase[from];			// R:: redundant?
	}

#if DO_MAXMIX
	nt->maxMix = alignTotal;
#endif
	return nt;
}

/////////////////////////////////

void
	ntConnect(Net *nt, Link *link)
{
	Cursor from = link->from;
	Cursor to   = link->to;
	Branch *br;
	Cursor bb;
	Cursor *up, *dn;

	if (nt->connection[from][to])
		return;
	nt->connection[from][to] = YES;
	nt->nLinks++;
	nt->nChildren[from]++;
	nt->nParents[to]++;
	if (nt->nParents[to] == 2)
		nt->nMixed++;

	bb = nt->freeBr;
	assert( bb != -1 );

	br = &nt->br[bb];
	nt->freeBr = br->nxtBr;
	br->fr = from;
	br->to = to;

#if SORTLINKS
	for (up = &nt->Ups[to]; *up != -1; up = &nt->br[*up].nxtUp) {
		if (nt->br[*up].fr > from)
			break;
	}
	br->nxtUp = *up;
	*up = bb;

	for (dn = &nt->Dns[from]; *dn != -1; dn = &nt->br[*dn].nxtDn) {
		if (nt->br[*dn].to > to)
			break;
	}
	br->nxtDn = *dn;
	*dn = bb;
#else 
	br->nxtUp = nt->Ups[to];
	nt->Ups[to] = bb;
	br->nxtDn = nt->Dns[from];
	nt->Dns[from] = bb;
#endif
	br->nxtBr = -2;
	//assert( nt->inuse[from] );
	//assert( nt->inuse[to] );
}

void
	ntDisconnectAll(Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor from;

	for (from = 0; from < tx->nTotal; from++)
		ASET(nt->connection[from], 0, tx->nTotal);
	ASET(nt->nChildren, 0, tx->nTotal);
	ASET(nt->nParents, 0, tx->nTotal);
	nt->nLinks = 0;
	nt->nMixed = 0;

	for (from = 0; from < 2*tx->nExtant-1; from++)
		nt->br[from].nxtBr = from+1;
	nt->br[from].nxtBr = -1;
	nt->freeBr = 0;
	ASET(nt->Ups, -1, tx->nTotal);
	ASET(nt->Dns, -1, tx->nTotal);
}

void
	ntDisconnect(Net *nt, Link *link)
{
	Cursor from = link->from;
	Cursor to   = link->to;
	Branch *br;
	Cursor *up, *dn;

	if (!nt->connection[from][to])
		return;
	nt->connection[from][to] = NO;
	--nt->nLinks;
	--nt->nChildren[from];
	--nt->nParents[to];
	if (nt->nParents[to] == 1)
		--nt->nMixed;

	br = 0;
	for (up = &nt->Ups[to]; *up != -1; up = &br->nxtUp) {
		br = &nt->br[*up];
		assert( br->to == to );
		if (br->fr == from)
			break;
	}
	assert( *up != -1 );

	for (dn = &nt->Dns[from]; *dn != -1; dn = &br->nxtDn) {
		br = &nt->br[*dn];
		assert( br->fr == from );
		if (br->to == to)
			break;
	}
	assert( *dn != -1 );

	assert( *up == *dn );
	*up = br->nxtUp;
	*dn = br->nxtDn;
	br->nxtBr = nt->freeBr;
	nt->freeBr = br - nt->br;
}

//////////////////////////
//
//  Constraints
//
int
	ntConstraints(Net *nt, FILE *fcon)
{
	Cursor from = TXNOT;
	char buf[20];
	enum states { doFrom, doTo, doTime, } state;
	char *deMix = getenv("DEMIX");
	int mixStratum = (deMix) ? atoi(deMix) : 0;    // Sentinel: OK for ROOT to not be mixed.

	state = doFrom;
	while (fscanf(fcon, "%20s", buf) == 1) {
		Cursor taxon = txFind(nt->taxa, buf);
	
		if (state == doTime && *buf != '<') {
			int stratum = atoi(buf);
			if (stratum <= mixStratum && from != TXNOT)
				nt->banMixed[from] = YES;

			// deMix taxa starting with an underscore
			if (txName(nt->taxa, from)[0] == '_')
				nt->banMixed[from] = YES;
			continue;
		}

		if (taxon == TXNOT) {
			if (!strchr("<>", *buf)) {
				fprintf(stderr, "Taxon expected, got: %s\n", buf);
				return NO;
			}
			state = (state != doTo) ? doTo : doFrom;
			continue;
		}

		if (state == doFrom) {
			from = taxon;
			state = doTime;
			continue;
		}

		nt->noanc[from][taxon] = YES;
	}
	return YES;
}
				
////////////////////////////
//
//

//	nt->ancestors[from][anc] => transpose of following:
//	nt->descendents[to][des] => taxon(des) is a descendent of taxon(to)

void
	ntPropagate(Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor maxTax = nt->maxTax;
	Cursor ii, jj, kk, mm;
	Cursor i0;
	Flag **d = nt->descendents;
	int nBr = 2*tx->nExtant;

	mm = (tx->nTotal + ALIGN) & ~ALIGN;
#if 0
	ACPY(d[0], nt->connection[0], mm * mm);
	for (ii = 0; ii < maxTax; ii++)
		d[ii][ii] = 1;
#else
	ASET(d[0], 0, mm * mm);
	for (ii = 0; ii < nBr; ii++) {
		Branch *br = &nt->br[ii];
		if (br->nxtBr == -2)
			d[br->fr][br->to] = 1;
	}
	for (ii = 0; ii < maxTax; ii++)
		d[ii][ii] = 1;
#endif

	mm = (maxTax+ALIGN) / (ALIGN+1);

	// Some d[ii][] are always zero, so skip them
	for (i0 = 0; i0 < maxTax; i0++) {
		if (nt->nChildren[i0] > 0)
			break;
	}

	// No dependencies, so skip.
	if (i0 == maxTax)
		return;

	// Warshall's Algorithm to compute a transitive closure
	for (kk = 0; kk < maxTax; kk++) {
		for (ii = i0; ii < maxTax; ii++) {
			if (d[ii][kk]) {
				unsigned long *di = (unsigned long *) d[ii];
				unsigned long *dk = (unsigned long *) d[kk];
				jj = mm;
				do {
					*di++ |= *dk++;
				} while (--jj > 0);
			}
		}
	}
}

/* Progeny: List of descendents from node 'to': */

struct progeny {
	Cursor *descends;
	int nD;
};

struct progeny *
	ntProgeny(Net *nt, Cursor to)
{
	Taxa *tx = nt->taxa;
	Cursor des, nD = 0;
	static struct progeny pgy[1];

	if (!pgy->descends) {
		newmem(pgy->descends,tx->nTotal);
	}
	for (des = 0; des < nt->maxTax; des++) {
		if (nt->descendents[to][des] && nt->inuse[des])
			pgy->descends[nD++] = des;
	}
	pgy->nD = nD;

	return pgy;
}

int
	ntProgConstrained(Net *nt, Cursor from, struct progeny *pgy)
{
	Cursor anc, des, nD;
	Cursor *ds;

	nD = pgy->nD;
	ds = pgy->descends;

	// Check the ancestor/descendent tables.
	for (anc = 0; anc < nt->maxTax; anc++) {
		if (!nt->descendents[anc][from])
			continue;
		if (!nt->inuse[anc])
			continue;
		for (des = 0; des < nD; des++) {
			if (nt->noanc[anc][ds[des]])
				return YES;
		}
	}
	return NO;
}

int
	ntConstrained(Net *nt, Cursor from, Cursor to)
{
	struct progeny *pgy;

	pgy = ntProgeny(nt, to);
	return ntProgConstrained(nt, from, pgy);
}

// NB: Incremental constraint code also in hypPoly().
int
	ntHypConstrained(Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor anc;
	Cursor  up, dn;

#if 0 && DO_CLINK
	// All of these should have taken care of earlier.
	Cursor tt;
	for (tt = 0; tt < tx->nExtant; tt++) {
		Cursor org = txCorrecting(tx,tt);
		if (org != TXNOT && !nt->connection[org][tt]) {
			return YES;
		}
	}
#endif

	// If we have an outgroup, don't let it be mixed.
	anc = nt->outgroup;
	if (anc != -1) {
		while (nt->nParents[anc] > 0) {
			if (nt->nParents[anc] > 1)
				return YES;
			up = nt->Ups[anc];
			anc = nt->br[up].fr;
		}
	}

	// Forbid hypothetical documents that have more
	// mixed children than unmixed, non-fragmentary kids.
	for (anc = tx->nExtant; anc < nt->maxTax; anc++) {
		int unmixed, mixed;
		if (nt->nChildren[anc] == 0)
			continue;
		if (!nt->inuse[anc])
			continue;
		mixed = unmixed = 0;

		for (dn = nt->Dns[anc]; dn != -1; dn = nt->br[dn].nxtDn) {
			Cursor des = nt->br[dn].to;

			if (nt->nParents[des] > 1)
				mixed++;
			else if (!txFrag(tx,des))
				unmixed++;
		}
		if (nt->nParents[anc] == 0)
			unmixed--;
		if (nt->nParents[anc] > 1)
			mixed++;
		if (mixed > unmixed)
			return YES;
	}

	return NO;
}


#if DO_POLE
static Length
	ntPoleStretch(Net *nt, Cursor to)
{
	// R:: Consider do-while() loop, and init hi,lo to first
	Length hi = 0L, lo = ~0L;
	for (Cursor up = nt->Ups[to]; up != ERR; up = nt->br[up].nxtUp) {
		Cursor p = nt->br[up].fr;
		Length pole = nt->poles[p];
		if (pole > hi) hi = pole;
		if (pole < lo) lo = pole;
	}
	return hi-lo;
}

int
	ntPoleCheck(Net *nt, Cursor to)
{
	if (nt->nParents[to] < 2)
		return YES;
#if NO_TRIPS
	if (nt->nParents[to] > 2)
		return NO;
#endif

	if (nt->banMixed[to])
		return NO;

	Length toPole = nt->poles[to];
	Length stretch = ntPoleStretch(nt, to);
	return stretch <= toPole;
}

int
	ntKidsPoleCheck(Net *nt, Cursor par)
{
	// Pole check all of par's kids now
	for (Cursor dn = nt->Dns[par]; dn != ERR; dn = nt->br[dn].nxtDn) {
		Cursor to = nt->br[dn].to;
		if (!ntPoleCheck(nt, to))
			return NO;
	}
	return YES;
}

void
	ntPolesOK(Net *nt)
{
	for (Cursor tt = 0; tt < nt->maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;
		Length pole = txPole(nt->taxa, tt);
		if (nt->poles[tt] != pole)
			printf("\nBad pole for tt:%d => txPole(tx,tt):"CST_F" != nt->poles[tt]:"CST_F"\n", tt, pole, nt->poles[tt]);
		assert( nt->poles[tt] == pole );
	}
}
#endif

///////////////
//
// Links helper
//

static Link *
	ntPreorder(Net *nt, Link *end, Cursor node)
{
	Taxa *tx = nt->taxa;
	Cursor dn;
	
	if (txVisit(tx,node))
		return end;
	txVisit(tx,node) = YES;

	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn) {
		end->from = node;
		end->to   = nt->br[dn].to;
		end++;
	}

	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn)
		end = ntPreorder(nt, end, nt->br[dn].to);

	return end;
}

static Cursor *
	ntPostorder(Net *nt, Cursor *end, Cursor node)
{
	Taxa *tx = nt->taxa;
	Cursor dn;
	
	if (txVisit(tx,node))
		return end;
	txVisit(tx,node) = YES;

	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn)
		end = ntPostorder(nt, end, nt->br[dn].to);
	*end++ = node;
	return end;
}

Cursor *
	ntPostorderNodes(Net *nt, Cursor **work)
{
	Taxa *tx = nt->taxa;
	Cursor tt;
	Cursor *end;

	if (!*work) {
		newmem(*work, tx->nTotal);
	}
	
	txUnvisit(tx);
	end = *work;
	for (tt = 0; tt < nt->maxTax; tt++) {
		if (!txVisit(tx,tt) && nt->inuse[tt])
			end = ntPostorder(nt, end, tt);
	}
	return end;
}

Link *
	ntPreorderLinks(Net *nt, Link **work)
{
	Taxa *tx = nt->taxa;
	Cursor tt;
	Link *end;

	if (!*work) {
		newmem(*work, nt->nLinks);
	}
	
	txUnvisit(tx);
	end = *work;
	for (tt = 0; tt < nt->maxTax; tt++) {
		if (nt->nParents[tt] == 0 && nt->inuse[tt])
			end = ntPreorder(nt, end, tt);
	}
	return end;
}

void
	ntRestoreLinks(Net *nt, Link *work, Link *end)
{
	Link *link;

	ntDisconnectAll(nt);
	for (link = work; link < end; link++)
		ntConnect(nt, link);
}

Cursor
	ntRoot(Net *nt)
{
	Cursor tt;

	for (tt = 0; tt < nt->maxTax; tt++) {
		if (nt->nParents[tt] == 0 && nt->inuse[tt])
			return tt;
	}
	fprintf(stderr, "Cannot find root, aborting...\n");
	abort();
	return TXNOT;
}
/////////////////////////////////////////////////////

/////////////////////////////////
//
// Cost Calculation
//

unsigned long TaxonCode;

struct memo {
	unsigned long fc;	// From code
	unsigned long tc;	// To code
	Length cume;		// Cumulative cost
	Flag *vc;			// Variant costs
};
enum Memoed { NA, MISS, HIT };

#define vcAlias(nt,i,vc) ((nt)->vcs[i]  = (vc))
#define vcUnalias(nt,i)  ((nt)->vcs[i]  = (nt)->vcBase[i])
#define vcClean(nt,i)    ((nt)->vcs[i] == (nt)->vcBase[i])

static void
	ntResetCost(Net *nt, struct memo **Memo)
{
	register Taxa *tx = nt->taxa;
	Cursor tt;
	struct memo *mm;

	if (*Memo)
		return;

	newmem(*Memo, tx->nTotal);
	mm = *Memo;
	for (tt = 0; tt < tx->nTotal; tt++) {
		vcUnalias(nt,tt);
		mm->fc = mm->tc = -1;
		mm->cume = -1;
		newmem(mm->vc, (tx->nVunits+ALIGN) &~ ALIGN);
		ASET(mm->vc, 0, (tx->nVunits+ALIGN) &~ ALIGN);
		mm++;
	}
}

Length
	ntCost(Net *nt)
{
	Cursor tt, root = TXNOT;
	register Taxa *tx = nt->taxa;
	Length cost = 0;
	Cursor maxTax = nt->maxTax;
	int ii, nBr = 2*tx->nExtant;
	struct memo *mm;
	static struct memo *Memo;

	ntResetCost(nt, &Memo);

	cost = 0;

	for (tt = 0; tt < maxTax; tt++) {
		if (!nt->inuse[tt])
			continue;

		mm = &Memo[tt];
		nt->cumes[tt] = tx->nVunits;
		cost += tx->nVunits;

		if (nt->nParents[tt] == 0) {
			// Assign Root to the first orphan
			if (root == TXNOT)
				root = tt;
			else {
				vcUnalias(nt,tt);
				memset(nt->vcs[tt], 1, tx->nVunits);
			}
		} else if (nt->nParents[tt] > 1) {
			Length ret_cost = RetCost * (nt->nParents[tt]-1);
			nt->cumes[tt] += ret_cost;
			cost += ret_cost;
		}
	}

	// Clear root's costs
	vcUnalias(nt,root);
	memset(nt->vcs[root], 0, tx->nVunits);
	cost -= nt->cumes[root];
	nt->cumes[root] = 0;

	// Loop over each link
	txUnvisit(tx);
	for (ii = 0; ii < nBr; ii++) {
		Branch *br = &nt->br[ii];
		Cursor fr, to;
		Length vcost = 0;
		enum Memoed memo;
		Length *cume;				// R:: Is this really faster using nt->cumes[to] ?
		int nPars;

		if (br->nxtBr != -2)
			continue;

		to = br->to;
		fr = br->fr;
		mm = &Memo[to];
		cume = &nt->cumes[to];

		nPars = nt->nParents[to];
		if (nPars != 1)
			memo = NA;
		else if (mm->tc == nt->codes[to] && mm->fc == nt->codes[fr])
			memo = HIT;
		else
			memo = MISS;

		if (memo == HIT) {
			assert( !txVisit(tx,to) );
			txVisit(tx,to) = YES;
			vcost = *cume - mm->cume;
			vcAlias(nt,to,mm->vc);
		} else {
			Flag *vc = vcUnalias(nt,to);

			// Loop over each character
			if (!txVisit(tx,to)) {
				vcost += stLinkCost(tx, fr, to, vc);	// May change *cume ?
				txVisit(tx,to) = YES;
			} else
				vcost += stRetLinkCost(tx, fr, to, vc);
		}
		assert( *cume - vcost >= 0 );
		assert( *cume - vcost <= tx->nVunits + nPars*RetCost );
		*cume -= vcost;
		cost -= vcost;

		// Memoize values
		if (memo == MISS) {
			mm->tc = nt->codes[to];
			mm->fc = nt->codes[fr];
			mm->cume = *cume;
			assert( vcClean(nt,to) );
			ACPY(mm->vc, nt->vcs[to], tx->nVunits);
		}
	}

	return cost;
}

/////////////////////////////////
//
// Incremental Cost Calculation
//

static Length
	ntIncSaveLink(Net *nt, Cursor from, Cursor to, struct progeny *pgy,
		Length cost, Length bound, Link *add)
{
	if (cost < bound
#if SORTLINKS
	|| (cost == bound && (from < add->from
	|| (from == add->from && to <= add->to)))
#endif
	) {
		Link link[1];

		// Slow check, do as late as possible
		if (ntProgConstrained(nt, from, pgy))
			return bound;
		link->from = from;
		link->to = to;
		ntConnect(nt, link);
		if (!ntHypConstrained(nt)) {
			*add = *link;
			bound = cost;
		}
		ntDisconnect(nt, link);
	}
	return bound;
}

// Inner loop for linking to a node with no parents
//		Assumes ntCost() was called upstream to set up basecost,
//		nt->cost[]s
static Length
	ntIncNewLink(Net *nt, Taxa *tx, Length basecost, Length bound,
		Cursor to, Link *add)
{
	Cursor from;
	unsigned vv;
	vunit *vTo, *vFrom;
	Length hiBound;
	struct progeny *pgy;

	pgy = ntProgeny(nt, to);
	vTo = txRdgs(tx,to);
	for (from = 0; from < nt->maxTax; from++) {
		Length cost;
		unsigned cm;	// Cost margin = cost - hiBound

		if (nt->noanc[from][to])
			continue;				// Quickie constraint check
		if (nt->connection[from][to])
			continue;				// Already linked.
		if (!nt->inuse[from])
			continue;

		cost = basecost;

		// Try to determine savings at the to-node.
		// ...if the max possible savings can get us below
		// ...the bound, then skip.
		hiBound = bound + tx->nVunits;
		if (cost >= hiBound)
			continue;

		// Critical loop
		cm = cost - hiBound;
		vv = tx->nVunits - 1;
		vFrom = txRdgs(tx,from);
		do {
			vunit toState = vTo[vv];
			vunit frState = vFrom[vv];
			cm += ((frState != toState) & (toState != MISSING));
			vv--;
		} while (vv < cm);
		cost = (int) cm + hiBound - tx->nVunits + 1 + (int) vv;

		/*	Some discussion of this critical loop is necessary.  The goal
		 *	is to use branch & bound without getting bogged the loop control.
		 *	Previously, vv < nV was tested separately from cost + vv >= hiB
		 * 	and, because these tests are for counters, they are merged into
		 *	a single unsigned comparison where the usually positive lhs
		 *	counts down to -1 and the usually negative rhs counts up to 0.
		 * 	NB: Some quirks of unsigned arithmetic are exploited.  The lhs
		 *	stops when vv > tx->nVunits-1, becoming -1 and forcing a FALSE.
		 *	The rhs has to be TRUE when	cost is clearly under the bound
		 *	but forces a FALSE when cost+vv == hiBound (since we're adding
		 *	a link, ties will always lose in the tie-breaker).
		 *
		 *	Also, the test has a strength reduction of nV - vv to vv, and
		 *	cost + vv - hiBound to cm.
		 *****************************************************************/

		bound = ntIncSaveLink(nt, from, to, pgy, cost, bound, add);
	}
	return bound;
}

// Inner loop for linking to a node with some parents
struct vdiff {
	Cursor vv;		// Vunit were there is a cost
	vunit toState;	// State of to-item
};

struct vdcache {
	CodeID nodeCode;
	CodeID parCode;
	struct vdiff *vdiffs;
	struct vdiff *vdend;
};

static struct vdiff *
	ntDiffLinks(Net *nt, Taxa *tx, Cursor to, struct vdiff **pvdend)
{
	static struct vdcache *vds;
	struct vdcache *vdc;
	struct vdiff *vdiffs;
	struct vdiff *vd;
	Cursor tt, vv;
	Flag *vc;
	vunit *vTo;

	if (!vds) {
		newmem(vds, tx->nTotal);
		for (tt = 0; tt < tx->nTotal; tt++) {
			vds[tt].nodeCode = vds[tt].parCode = -1;
			vds[tt].vdiffs = vds[tt].vdend = 0;
		}
	}

	vdc = &vds[to];
	if (!vdc->vdiffs) {
		newmem(vdc->vdiffs, tx->nVunits);
	}
	vdiffs = vdc->vdiffs;

	// Check vdiff cache
	if (nt->nParents[to] == 1) {
		Cursor up = nt->Ups[to];
		Cursor fr = nt->br[up].fr;
		if (vdc->nodeCode == nt->codes[to] && vdc->parCode == nt->codes[fr]) {
			*pvdend = vdc->vdend;
			return vdiffs;	
		}
		vdc->nodeCode = nt->codes[to];
		vdc->parCode = nt->codes[fr];
	} else
		vdc->nodeCode = vdc->parCode = -1;

	// Find the relative few vunits where there is a cost
	// Assumes ntCost() was called upstream to set up basecost, nt->cost[]s
	vd = vdiffs;
	vc = nt->vcs[to];
	vTo = txRdgs(tx,to);
	for (vv = 0; vv < tx->nVunits; vv += ALIGN+1) {
		if (*(unsigned long *)&vc[vv] == 0)
			continue;

#define SET_VDIFF(vc, vd, _vv, vt) \
	if (vc[_vv]) { \
		vd->vv = (_vv); \
		vd->toState = vt[_vv]; \
		vd++; \
	} else

		SET_VDIFF(vc, vd, vv+0, vTo);
		SET_VDIFF(vc, vd, vv+1, vTo);
		SET_VDIFF(vc, vd, vv+2, vTo);
		SET_VDIFF(vc, vd, vv+3, vTo);
#if defined (_LP64)
		// I want #if ALIGN == 7, but ALIGN contains sizeof and cannot be used here.
		// So, I'll check the architecture and sanity check with a "static assert" to make sure.
switch (0) case 0: case (ALIGN - 3): ;	// Make sure the word size is not actually 4.
		SET_VDIFF(vc, vd, vv+4, vTo);
		SET_VDIFF(vc, vd, vv+5, vTo);
		SET_VDIFF(vc, vd, vv+6, vTo);
		SET_VDIFF(vc, vd, vv+7, vTo);
#endif

	}
	*pvdend = vd;

	// Cache the end.
	vdc->vdend = vd;
	return vdiffs;
}

static Length
	ntIncRetLink(Net *nt, Taxa *tx, Length basecost, Length bound,
		Cursor to, Link *add)
{
	Cursor from;
	vunit *vFrom;
	struct vdiff *vdiffs;
	struct vdiff *vd, *vdend;
	struct progeny *pgy;
	int nV;
#if DO_POLE
	Length hi = 0L, lo = ~0L;
	Length toPos = 0L;
	Length frPos = 0L;

#if NO_TRIPS
	if (nt->nParents[to] > 1)
		return bound;
#endif
	if (nt->banMixed[to])
		return bound;
#endif

	pgy = ntProgeny(nt, to);

	vdiffs = ntDiffLinks(nt, tx, to, &vdend);
	nV = vdend - vdiffs;

#if DO_POLE
	for (Cursor up = nt->Ups[to]; up != ERR; up = nt->br[up].nxtUp) {
		Cursor p = nt->br[up].fr;
		Length pole = nt->poles[p];
		if (pole > hi) hi = pole;
		if (pole < lo) lo = pole;
	}
	toPos = nt->poles[to];
#endif

	for (from = 0; from < nt->maxTax; from++) {
		Length cost, hiBound;
		unsigned int cm, vv;

		if (nt->noanc[from][to])
			continue;				// Quickie constraint check
		if (nt->connection[from][to])
			continue;				// Already linked.
		if (!nt->inuse[from])
			continue;

#if DO_POLE
		frPos = nt->poles[from];
		Length hi2 = (frPos > hi) ? frPos : hi;
		Length lo2 = (frPos < lo) ? frPos : lo;
		if (hi2 - lo2 > toPos)
			continue;
#endif
		cost = basecost;
		// Try to determine savings at the to-node.
		// ...if the max possible savings can get us below
		// ...the bound, then skip.
		cost += RetCost;
		hiBound = bound + nV;
		if (cost > hiBound)
			continue;

		// Critical loop
		vFrom = txRdgs(tx,from);

		cm = cost - hiBound;
		vv = nV - 1;
		if (vv != -1) do {
			vd = vdiffs + vv;
			vv--;
			if (vFrom[vd->vv] != vd->toState)
				cm++;
		} while (vv < cm);
		cost = (int) cm + hiBound - (nV - (int) vv - 1);

		bound = ntIncSaveLink(nt, from, to, pgy, cost, bound, add);
	}
	return bound;
}

			// Add a link to the Net
Length
	ntIncLink(Net *nt, Length basecost, Length bound, Cursor dst, Link *add)
{
	Taxa *tx = nt->taxa;

	assert( dst != TXNOT );

	if (nt->nParents[dst] == 0)
		bound = ntIncNewLink(nt, tx, basecost, bound, dst, add);
	else
		bound = ntIncRetLink(nt, tx, basecost, bound, dst, add);

	return bound;
}
