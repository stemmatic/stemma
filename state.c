
/*
	state.c - routines for calculating character state
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "state.h"

typedef unsigned char States;

typedef vunit *Kin;

int
	stKids(Net *nt, Cursor node, vunit **rdgs, Kin **pKids)
{
	Taxa *tx = nt->taxa;
	Cursor dn;
	Kin *kids = *pKids;

	if (!kids) {
		newmem(kids, tx->nTotal*2);
		*pKids = kids;
	}

	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn) {
		Cursor kk = nt->br[dn].to;
		*kids++ = rdgs[kk];
	}
	*kids = 0;
	return kids - *pKids;
}

int
	stFolks(Net *nt, Cursor node, vunit **rdgs, Kin **pFolks)
{
	Taxa *tx = nt->taxa;
	Cursor up;
	Kin *folks = *pFolks;

	if (!folks) {
		newmem(folks, tx->nTotal*2);
		*pFolks = folks;
	}

	for (up = nt->Ups[node]; up != -1; up = nt->br[up].nxtUp) {
		Cursor t = nt->br[up].fr;
		*folks++ = rdgs[t];
	}
	*folks = 0;
	return folks - *pFolks;
}

// Returns nInLaws

int
	stInLaws(Net *nt, Cursor node, vunit **rdgs, Kin **pInLaws)
{
	Taxa *tx = nt->taxa;
	Cursor dn;
	Kin *inLaws = *pInLaws;
	int nInLaws = 0;

	if (!inLaws) {
		newmem(inLaws, tx->nTotal*2);
		*pInLaws = inLaws;
	}

	for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn) {
		Cursor kk = nt->br[dn].to;

		*inLaws++ = rdgs[kk];

		// Is the value for taxon kk spoken for by another parent?
		if (nt->nParents[kk] > 1) {
			Cursor up;
			for (up = nt->Ups[kk]; up != -1; up = nt->br[up].nxtUp) {
				Cursor t = nt->br[up].fr;
				if (t != node)
					*inLaws++ = rdgs[t];
			}
			nInLaws++;
		}
		*inLaws++ = 0;
	}
	*inLaws = 0;
	return nInLaws;
}

// Compute state of hypothetical Ancestor
//
// A -> B -> C,D  ==> A -> [1] -> B,C,D
//
// Determine the ancestral state by scoring the states of
// A, C, D, and picking the best state, breaking ties in
// favor of the parent.  (Except: all MISSING kids -> MISSING)
//	

/*	Helper routines for stFixNode():
	stSetRet() - Most General Purpose code and slowest, used for mixed kids
	stSetSimp() - Semi-special case: Simple Kids, simple parent
*/
static void
	stSetRet(Net *nt, Cursor node, Kin *folks, Kin *inLaws)
{
	Taxa *tx = nt->taxa;
	vunit *base = txBase(tx,node);
	int nPars = nt->nParents[node];
	int vv;

	for (vv = 0; vv < tx->nVunits; vv++) {
		static States score[MAXSTATES];
		int total = nt->nChildren[node];
		unsigned int isPar;
		unsigned int maxS;
		unsigned int parS;
		unsigned int sBest = MISSING;
		Kin  *fol, *inl;
		vunit *pb, *kb;
		int nMiss = 0;

		ZERO(score, dimof(score));
		isPar = 0;
		for (fol = folks; (pb = *fol); fol++) {
			vunit pv = pb[vv];
			score[pv] = 1;
			isPar |= (1 << pv);
			if (pv-1 < sBest-1)
				sBest = pv;
		}
		maxS = (nPars > 0);
		parS = maxS;

		for (inl = inLaws; (kb = *inl); inl++) {
			vunit *tb;
			vunit kv = kb[vv];
			
			for (inl++; (tb = *inl); inl++) {
				if (tb[vv] == kv)
					kv = MISSING;
			}
			if (kv != MISSING) {
				States sc = ++score[kv];
				int parK = (isPar & (1 << kv)) != 0;

				if (sc > maxS) {
					maxS = sc;
					parS = parK;
					sBest = kv;
				} else if (sc < maxS)
					;
				else if (parK > parS) {
					parS = parK;
					sBest = kv;
				} else if (parK < parS)
					;
				else if (kv-1 < sBest-1)
					sBest = kv;
			} else
				nMiss++;
		}
		if (nMiss == total)
			sBest = MISSING;

		base[vv] = sBest;
	}
}

/*
	stSetSimp - Semi-general purpose code for setting the state of
	simple hypotheticals.  No reticulation on node or kids.

*/

static void
	stSetSimp(Net *nt, Cursor node, vunit *pb, Kin *vkids)
{
	Taxa *tx = nt->taxa;
	vunit *base = txBase(tx,node);
	int nKids = nt->nChildren[node];
	int vv;

	for (vv = 0; vv < tx->nVunits; vv++) {
		static States score[MAXSTATES];
		Kin *kids;
		vunit *kb;
		vunit kv, pv = pb[vv];
		unsigned int isPar = (1 << pv);
		unsigned int maxS = 1, parS = 1;
		unsigned int sBest = pv;
		int nMiss = 0;

		ZERO(score, dimof(score));
		score[pv] = 1;

		for (kids = vkids; (kb = *kids); kids++) {
			kv = kb[vv];
			
			if (kv != MISSING) {
				States sc = ++score[kv];
				int parK = (isPar & (1 << kv)) != 0;

				if (sc > maxS) {
					maxS = sc;
					parS = parK;
					sBest = kv;
				} else if (sc < maxS)
					;
				else if (parK > parS) {
					parS = parK;
					sBest = kv;
				} else if (parK < parS)
					;
				else if (kv-1 < sBest-1)
					sBest = kv;
			} else
				nMiss++;
		}
		if (nMiss == nKids)
			sBest = MISSING;

		base[vv] = sBest;
	}
}

/////////////////////////////////////////////////////////
// stFixNode - recalculate states for an entire node

void
	stFixNode(Net *nt, Cursor node)
{
	Taxa *tx = nt->taxa;
	static Kin *folks = 0;
	static Kin *vkids = 0;
	static Kin *inLaws = 0;
	int nKids, nPars;
	int isSimple;

	if (node < tx->nExtant)
		return;

	nKids = nt->nChildren[node];
	isSimple = (stInLaws(nt, node, &txBase(tx,0), &inLaws) == 0);
	nPars = stFolks(nt, node, &txBase(tx,0), &folks);
	assert( nPars == nt->nParents[node] );

	if (nKids == 0) {
		int vv;
		for (vv = 0; vv < tx->nVunits; vv++)
			txBase(tx,node)[vv] = MISSING;
	} else if (isSimple && nPars == 1) {
		vunit *bp = folks[0];
		vunit *bn = txBase(tx,node);
		int vv;

		stKids(nt, node, &txBase(tx,0), &vkids);

		if (nKids == 3) {
			vunit *b1 = vkids[0];
			vunit *b2 = vkids[1];
			vunit *b3 = vkids[2];

			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit s1 = b1[vv];
				vunit s2 = b2[vv];
				vunit s3 = b3[vv];
				vunit sp = bp[vv];
				bn[vv] = (s1 == s2 && s2 == s3) ? s1
					: (sp != MISSING && (sp == s1 || sp == s2 || sp == s3)) ? sp
					: (s1 != MISSING && (s1 == s2 || s1 == s3)) ? s1
					: (s2 != MISSING && (s2 == s3)) ? s2 : sp;
			}
		} else if (nKids == 2) {
			vunit *b1 = vkids[0];
			vunit *b2 = vkids[1];

			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit s1 = b1[vv];
				vunit s2 = b2[vv];
				vunit sp = bp[vv];
				bn[vv] = (s1 == s2) ? s1 : sp;
			}
		} else if (nKids == 1) {
			vunit *b1 = vkids[0];

			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit s1 = b1[vv];
				vunit sp = bp[vv];
				bn[vv] = (s1 == MISSING) ? s1 : sp;
			}
		} else
			stSetSimp(nt, node, folks[0], vkids);
	} else // call general case
		stSetRet(nt, node, folks, inLaws);

	txPermuteTaxon(tx, node);
	nt->codes[node] = TaxonCode++;
}

/////////////////////////////////////////////////////////////
//
//	Incremental Cost Calculation integrated with State Calculation

static Length
	stAncGenExtantCost(Net *nt, Cursor hyp, Kin *folks, Kin *inLaws,
		Length cost, Length bound)
{
	const Taxa *tx = nt->taxa;
	const int nKids = nt->nChildren[hyp];
	const int nPars = nt->nParents[hyp];
	const int isRoot = (nPars == 0), oneP = (nPars == 1);
	vunit *hb;
	int vv;

	/* Calculate cost for extant (non-hypothetical) node */
	hb = txRdgs(tx,hyp);

	// Extant, no Kids
	if (nKids == 0) {

		// No Kids, no Parents
		if (isRoot)
			return cost;

		// No Kids, one parent is done upstairs
		assert( !oneP );

		// No Kids, at least two parents
		for (vv = 0; vv < tx->nVunits; vv++) {
			Kin *fol;
			vunit *pb;
			vunit hv;
			int pc = YES;

			hv = hb[vv];
			if (hv == MISSING)
				continue;
			fol = folks;
			for (fol = folks; (pb = *fol); fol++)
				pc &= (pb[vv] != hv);
			if (pc)
				cost++;
		}
		return cost;
	}

	// Extant Calc Code with one Parent
	if (oneP) {
		vunit *pb = folks[0];
		for (vv = 0; vv < tx->nVunits; vv++) {
			Kin *inl;
			vunit *kb;
			vunit hv;

			hv = hb[vv];
			for (inl = inLaws; (kb = *inl); inl++) {
				vunit *tb;
				vunit kv = kb[vv];
				int kc = (hv != kv) & (kv != MISSING);

				for (inl++; (tb = *inl); inl++)
					kc &= (tb[vv] != kv);

				if (kc)
					cost++;
			}
			if (pb[vv] != hv && hv != MISSING)
				cost++;   // assert( nPars == 1 );
		}
		return cost;
	}

	// General Purpose Extant Calc Code
	for (vv = 0; vv < tx->nVunits; vv++) {
		Kin *inl;
		Kin *fol;
		vunit *kb, *pb;
		vunit hv;
		int pc = YES;

		hv = hb[vv];
		for (inl = inLaws; (kb = *inl); inl++) {
			vunit *tb;
			vunit kv = kb[vv];
			int kc = (hv != kv) & (kv != MISSING);

			for (inl++; (tb = *inl); inl++)
				kc &= (tb[vv] != kv);

			if (kc) {
				cost++;
				if (cost > bound)
					return cost;
			}
		}
		if (hv == MISSING)
			continue;
		for (fol = folks; (pb = *fol); fol++)
			pc &= (pb[vv] != hv);
		if ((fol != folks) && pc)
			cost++;
	}
	return cost;
}

static Length
	stAncGenRootCost(Net *nt, Cursor hyp, Kin *vkids, Kin *inLaws,
		int isSimple, Length cost)
{
	const Taxa *tx = nt->taxa;
	const int nKids = nt->nChildren[hyp];
	static States score[MAXSTATES];
	States total;
	int vv;

	if (isSimple && nKids == 1) {
		vunit *kb1 = vkids[0];

		for (vv = 0; vv < tx->nVunits; vv++)
			cost += (kb1[vv] != MISSING);
		return cost;
	}

	// Hyp, isRoot, mixed kids or 3+ kids
	for (vv = 0; vv < tx->nVunits; vv++) {
		Kin *inl;
		vunit *kb;
		States maxS = 0;

		ZERO(score, dimof(score));
		total = 1;					// Get kid contribs later
		
		for (inl = inLaws; (kb = *inl); inl++) {
			vunit *tb;
			vunit kv = kb[vv];
			
			for (inl++; (tb = *inl); inl++) {
				if (tb[vv] == kv)
					kv = MISSING;
			}
			total++;						// Get kid's contrib
			score[kv]++;
			if (kv != MISSING && score[kv] > maxS)
				maxS = score[kv];
		}
		total -= score[MISSING];
		cost += total - maxS;
	}
	cost -= tx->nVunits;

	return cost;
}

static Length
	stAncGenCost(Net *nt, Cursor hyp, Kin *inLaws, Kin *folks, Length cost)
{
	const Taxa *tx = nt->taxa;
	const int nPars = nt->nParents[hyp];
	static States score[MAXSTATES];
	States total;
	int vv;

	// General Purpose Hypothetical Node Code:
	for (vv = 0; vv < tx->nVunits; vv++) {
		Kin *fol, *inl;
		vunit *kb, *pb;
		States maxS = (nPars > 0);

		ZERO(score, dimof(score));
		for (fol = folks; (pb = *fol); fol++) {
			vunit pv = pb[vv];
			score[pv] = 1;
		}
		total = 1;
		score[MISSING] = 0;
		
		for (inl = inLaws; (kb = *inl); inl++) {
			vunit *tb;
			vunit kv = kb[vv];
			
			for (inl++; (tb = *inl); inl++) {
				if (tb[vv] == kv)
					kv = MISSING;
			}
			total++;
			score[kv]++;
			if (kv != MISSING && score[kv] > maxS)
				maxS = score[kv];
		}
		total -= score[MISSING];
		cost += total - maxS;
	}
	if (nPars == 0)
		cost -= tx->nVunits;

	return cost;
}

Length
	stAncCost(Net *nt, Cursor hyp, Length *cumes, Length cost, Length bound)
{
	const Taxa *tx = nt->taxa;
	int dn, vv;
	static Kin *inLaws = 0;
	static Kin *folks = 0;
	static Kin *vkids = 0;
	static States score[MAXSTATES];
	States total;
	vunit *pb;

	// Compile local topology into more easily accessible data structures
	int nKids = stKids(nt, hyp, &txRdgs(tx,0), &vkids);
	int nPars = stFolks(nt, hyp, &txRdgs(tx,0), &folks);
	int isRoot = (nPars == 0), oneP = (nPars == 1);
	int isSimple = (stInLaws(nt, hyp, &txRdgs(tx,0), &inLaws) == 0);

	assert( nKids == nt->nChildren[hyp] );
	assert( nPars == nt->nParents[hyp] );

	// Recalc. kids
	for (dn = nt->Dns[hyp]; dn != -1; dn = nt->br[dn].nxtDn) {
		Cursor kk = nt->br[dn].to;
		Length vunit_costs = 0;
		vunit_costs = cumes[kk];
		if (nt->nParents[kk] > 1)
			vunit_costs -= RetCost * (nt->nParents[kk]-1);
		cost -= vunit_costs;
	}
	if (cost > bound)
		return cost;

	// Mop up the new hypothetical node costs.
	if (nPars > 1)
		cost += RetCost * (nPars-1);
	if (cost > bound)
		return cost;

	if (hyp < tx->nExtant) {
		// Most common case: 0 kids, 1 parent
		if (nKids == 0 && oneP) {
			vunit *pb = folks[0];
			vunit *hb = txRdgs(tx,hyp);
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit hv = hb[vv];
				if (pb[vv] != hv && hv != MISSING)
					cost++;
			}
			return cost;
		}

		return stAncGenExtantCost(nt, hyp, folks, inLaws, cost, bound);
	}

	// No parent cases
	if (isRoot) {
		if (isSimple && nKids == 2) {
			vunit *kb1 = vkids[0];
			vunit *kb2 = vkids[1];

			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit kv1 = kb1[vv];
				vunit kv2 = kb2[vv];
				if (kv1 != kv2 && kv1 != MISSING && kv2 != MISSING)
					cost++;
			}
			return cost;
		}

		return stAncGenRootCost(nt, hyp, vkids, inLaws, isSimple, cost);
	}

	if (!oneP)
		return stAncGenCost(nt, hyp, inLaws, folks, cost);

	// One Parent (about 80% of cases):
	pb = folks[0];

	// Hyp, One Parent, Mixed Kids 
	if (!isSimple) {
		for (vv = 0; vv < tx->nVunits; vv++) {
			Kin *inl;
			vunit *kb;
			vunit pv;
			States maxS = 1;

			ZERO(score, dimof(score));
			pv = pb[vv];
			if (pv != MISSING)
				score[pv] = nPars;		// assert( nPars==1 );
			total = nPars;				// Add kids' contrib later
			
			for (inl = inLaws; (kb = *inl); inl++) {
				vunit *tb;
				vunit kv = kb[vv];
				
				for (inl++; (tb = *inl); inl++) {
					if (tb[vv] == kv)
						kv = MISSING;
				}
				total++;			// Kid's contribution
				score[kv]++;
				if (kv != MISSING && score[kv] > maxS)
					maxS = score[kv];
			}
			total -= score[MISSING];
			cost += total - maxS;
		}
		return cost;
	}

	// Hyp, One Parent, Simple, 2 kids
	if (nKids == 3) {
		vunit *kb1 = vkids[0];
		vunit *kb2 = vkids[1];
		vunit *kb3 = vkids[2];

		for (vv = 0; vv < tx->nVunits; vv++) {
			vunit pv = pb[vv];
			vunit kv1 = kb1[vv];
			vunit kv2 = kb2[vv];
			vunit kv3 = kb3[vv];

			if (kv1 != pv && kv1 != MISSING)
				cost++;
			if (kv2 != pv && kv2 != kv1 && kv2 != MISSING)
				cost++;
			if (kv3 != MISSING && ((kv3 != pv)
				? ((kv3 != kv1 && kv3 != kv2)
					|| (kv3 == kv1 && kv2 == pv && pv != MISSING)
					|| (kv3 == kv2 && kv1 == pv && pv != MISSING))
				: (kv3 != kv1 && kv1 == kv2 && kv1 != MISSING)))
				cost++;
		}
		return cost;
	}

	// Hyp, One Parent, Simple w/4+ kids.
	if (nKids > 3) {
		for (vv = 0; vv < tx->nVunits; vv++) {
			Kin *kids;
			vunit *kb;
			vunit pv;
			States maxS = 1;
			int nMiss = 0;

			ZERO(score, dimof(score));
			pv = pb[vv];
			if (pv != MISSING)
				score[pv] = nPars;		// assert( nPars==1 );
			total = nPars;				// Add kids' contrib later
			
			for (kids = vkids; (kb = *kids); kids++) {
				vunit kv = kb[vv];
				
				total++;			// Kid's contribution
				if (kv == MISSING)
					nMiss++;
				else {
					score[kv]++;
					if (score[kv] > maxS)
						maxS = score[kv];
				}
			}
			cost += total - maxS - nMiss;
		}
		return cost;
	}

	// Hyp, One Parent, Simple, 2 kids
	if (nKids == 2) {
		vunit *kb1 = vkids[0];
		vunit *kb2 = vkids[1];

		for (vv = 0; vv < tx->nVunits; vv++) {
			vunit pv = pb[vv];
			vunit kv1 = kb1[vv];
			vunit kv2 = kb2[vv];

			if (kv1 != pv && kv1 != MISSING)
				cost++;
			if (kv2 != kv1 && kv2 != pv && kv2 != MISSING)
				cost++;
		}
		return cost;
	}

	// Hyp, One Parent, One Simple Kid
	if (nKids == 1) {
		vunit *kb = vkids[0];
		for (vv = 0; vv < tx->nVunits; vv++) {
			vunit kv = kb[vv];
			vunit pv = pb[vv];
			if (pv != kv && kv != MISSING)
				cost++;
		}
		return cost;
	}

	// Hyp, One Parent, Zero simple kids
	return cost;
}

static Length
	stPolyRetCost(Net *nt, Cursor hyp, Cursor tt, Length cost)
{
	const Taxa *tx = nt->taxa;
	// NB: I could have much cleaner code for one in-law case.

	static Kin *inLaws = 0;
	Kin *inl;
#if DO_MP2
	vunit a1, a2;		// Actual states for kids
	int nMiss1 = 0, nMiss2 = 0;
#endif
	int vv;

	stInLaws(nt, hyp, &txRdgs(tx,0), &inLaws);

	for (vv = 0; vv < tx->nVunits; vv++) {
		vunit sp = txRdgs(tx,tt)[vv];
		vunit s1, s2, sh;
		vunit *kb, *tb;
		inl = inLaws;

		// Do first kid
		kb = *inl;
		s1 = kb[vv];
#if DO_MP2
		a1 = s1;
#endif
		for (inl++; (tb = *inl); inl++) {
			if (tb[vv] == s1)
				s1 = MISSING;
		}
		inl++;

		// Do second kid
		kb = *inl;
		s2 = kb[vv];
#if DO_MP2
		a2 = s2;
#endif
		for (inl++; (tb = *inl); inl++) {
			if (tb[vv] == s2)
				s2 = MISSING;
		}

		sh = (s1 == MISSING) ? s2 : (s2 == MISSING) ? s1 :
			 (1 + (s1==sp) >= 1 + (s2==sp)) ? s1 : s2;
		if (sh != sp && sh != MISSING)
			cost++;
		if (s1 != sh && s1 != MISSING)
			cost++;
		if (s2 != sh && s2 != MISSING)
			cost++;
#if DO_MP2
		if (sh == MISSING) {
			if (a1 != MISSING)
				nMiss1++;
			if (a2 != MISSING)
				nMiss2++;
		}
#endif
	}
#if DO_MP2
	// Punt on MP2 constraint violation
	if (nMiss1 > MaxMP2 || nMiss2 > MaxMP2)
		return -1;
#endif
	return cost;
}

Length
	stPolyCost(Net *nt, Cursor hyp, Cursor tt, Cursor t1, Cursor t2,
		Length *cumes, Length cost, Length bound)
{
	const Taxa *tx = nt->taxa;
	int vv;
	int noret = NO;
	const int nPar1 = nt->nParents[t1];
	const int nPar2 = nt->nParents[t2];
	vunit *vp, *v1, *v2;

	// Compile local topology into more easily accessible data structures
	if (nPar1 == 1 && nPar2 == 1)
		noret = YES;

	cost -= cumes[t1];
	if (nPar1 > 1)
		cost += RetCost * (nt->nParents[t1]-1);
	cost -= cumes[t2];
	if (nPar2 > 1)
		cost += RetCost * (nt->nParents[t2]-1);
	if (cost > bound)
		return cost;

	if (!noret)
		return stPolyRetCost(nt, hyp, tt, cost);

	vp = txRdgs(tx,tt);
	v1 = txRdgs(tx,t1);
	v2 = txRdgs(tx,t2);
	for (vv = 0; vv < tx->nVunits; vv++) {
		vunit sp = vp[vv];
		vunit s1 = v1[vv];
		vunit s2 = v2[vv];
		
		if (s1 != sp && s1 != MISSING)
			cost++;
		if (s2 != s1 && s2 != sp && s2 != MISSING)
			cost++;
	}
	return cost;
}

#if ROOTFIX
/*
	stPolyRootCost - predict the rootCost for a given POLY operation.

	With w -> a, b, c;  If a,b,c == 1,1,0, then w=1, for a rootCost=1
	with w -> [h], c; [h]->a,c; the same states mean w=0, for a rootCost=0
*/
Length
	stPolyRootCost(Net *nt, Cursor poly, Cursor parent, Cursor root,
		vunit *vRoot, Length rootCost)
{
	const Taxa *tx = nt->taxa;
	static Kin *rootKin = 0;
	static Kin *polyKin = 0;
	vunit *vA, *vB, *vC;		// The three root kids, A&B joined
	Length newRC;
	int vv;

	if (parent != root)
		return rootCost;

	if (nt->outgroup != TXNOT)
		return rootCost;

	if (root < tx->nExtant)
		return rootCost;

	// Simplest case: Root has two simple kids, and so does hyp.
	if (nt->nChildren[root] != 2)
		return rootCost;
	if (stInLaws(nt, root, &txRdgs(tx,0), &rootKin) > 0)
		return rootCost;

	assert( nt->nChildren[poly] == 2 );		// because we're Poly
	if (stInLaws(nt, poly, &txRdgs(tx,0), &polyKin) > 0)
		return rootCost;
	
	// Identify the three kids
	vA = polyKin[0]; assert( vA != 0 );
	vB = polyKin[2]; assert( vB != 0 );
	
	vC = (rootKin[0] == txRdgs(tx,poly)) ? rootKin[2] :
	     (rootKin[2] == txRdgs(tx,poly)) ? rootKin[0] : 0;
	assert( vC != 0 );

	// Calculate the root cost
	newRC = 0;
	for (vv = 0; vv < tx->nVunits; vv++) {
		int sA = vA[vv], sB = vB[vv], sC = vC[vv], sR = vRoot[vv];
		if (sC != 1 && ((sA == sB) ? sA : sR) != 1)
			newRC++;
	}

	return newRC;
}
#endif

/* Node Cost Routines */

int
	stLinkCost(Taxa *tx, Cursor fr, Cursor to, Flag *vc)
{
	vunit *vFrom = txRdgs(tx,fr);
	vunit *vTo = txRdgs(tx,to);
	unsigned vv = tx->nVunits - 1;
	int vcost = 0;

	//  First time here, so init vc based on vFrom
	//	For each link, zero step of vTo if any vFrom has same
	//	state.  Missing vTo characters are always free.
	//  For ?->0 transitions, set tc to -1, which sets vc[] to 2.
	do {
		vunit tv = vTo[vv];
		vunit fv = vFrom[vv];
		int tMiss = (tv == MISSING);
		int tc = ((fv == tv) | tMiss);

		vcost += tc;
		vc[vv] = 1 - tc;
	} while (vv-- > 0);
	return vcost;
}

int
	stRetLinkCost(Taxa *tx, Cursor fr, Cursor to, Flag *vc)
{
	vunit *vFrom = txRdgs(tx,fr);
	vunit *vTo = txRdgs(tx,to);
	unsigned vv = tx->nVunits - 1;
	int vcost = 0;

	do {
		vunit tv = vTo[vv];
		Flag vcv = vc[vv];

		// Second time through, all the to-MISSINGs have vcv==0
		if (vcv == 0)
			/* Variant already accounted for */;
		else if (vFrom[vv] == tv) {
			vcost++;
			vc[vv] = 0;
		}
	} while (vv-- > 0);
	return vcost;
}
