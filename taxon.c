/*
	taxon.c - Routines for handling taxa
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "stemma.h"
#include "taxon.h"

// Prefix of taxon names to ignore the singularity count.
#define IGNORE_COUNT_PREFIX '_'

static vunit *
	txInVars(Taxa *tx, Cursor tt, const char *vars)
{
	vunit *v, *start;
	int nVunits = tx->nVunits;
	int nMissing = 0;
	char *st0 = STATES;

	v = start = tx->base[tt];

	while (v - start < nVunits) {
		vunit var = (vars) ? *vars++ : '?';
		char *st1 = strchr(st0, var);

		if (!st1) {
			fprintf(stderr, "Invalid state '%c' for taxon %s at vunit %d.\n",
				var, tx->taxa[tt].name, (int) (v - start));
			assert( st1 != 0 );
		}
		if (var == '?')
			nMissing++;
		*v++ = st1 - st0;
	}

	txFrag(tx,tt) = (vars) ? (nMissing > nVunits/3) : NO;

	return start;
}

static enum VarType
	txVarType(Taxa *tx, Cursor vv)
{
	static vunit states[MAXSTATES];	// States attested among kids
	static int count[MAXSTATES];	// Count of corresponding state
	Cursor ss, nStates=0;			// Number of different states
	int   dblCount=0;				// Count for the twice attested
	vunit majState = '?';
	int	  majCount = 0;			    // Majority state, count

	Cursor tt;						// Taxa: taxon

	// Collect state counts
	for (tt = 0; tt < tx->nExtant; tt++) {

		// Skip missing states
		if (txRdgs(tx,tt)[vv] == MISSING)
			continue;
		
		// Find state index, adjust count if found, else add one
		for (ss = 0; ss < nStates; ss++) {
			if (states[ss] == txRdgs(tx,tt)[vv])
				break;
		}
		if (ss == nStates) {
			ss = nStates++;
			states[ss] = txRdgs(tx,tt)[vv];
			if (*txName(tx,tt) != IGNORE_COUNT_PREFIX)
				count[ss] = 1;
		} else if (*txName(tx,tt) != IGNORE_COUNT_PREFIX) {
			count[ss]++;
			if (count[ss] == 2)
				dblCount++;
		}
	}

	assert( nStates <= MAXSTATES );
	if (nStates <= 1)
		return Constant;

	// Upweight singulars
	for (ss = 0; ss < nStates; ss++) {
		if (count[ss] > majCount) {
			majCount = count[ss];
			majState = states[ss];
		}
		if (count[ss] == 1) {
			for (tt = 0; tt < tx->nExtant; tt++) {
				if (txRdgs(tx,tt)[vv] == states[ss]) {
					txIsSingular(tx,tt)[vv] = YES;
					txNSings(tx,tt)++;
				}
			}
		}
	}

	tx->majority[vv] = majState;
	if (dblCount <= 1)
		return Singular;

	return Informative;
}

Taxa *
	txScan(FILE *fp)
{
	Cursor nt, vv;
	static char name[TAXNAME];
	char *vbuf;
	Taxa *taxa;
	int alignTotal, alignVunits;

	newmem(taxa, 1);

	// Scan number of taxa and variation units.
	if (fscanf(fp, "%d %d", &taxa->nExtant, &taxa->nVunits) != 2)
		return 0;

		// +1 for temporary hypothetical if there are no reticulations
	taxa->nTotal = taxa->nExtant + (taxa->nExtant-1) + 1;

	// Compute aligned for fostering 32-bit memory moves
	alignTotal  = (taxa->nTotal +ALIGN) & ~ALIGN;
	alignVunits = (taxa->nVunits+ALIGN) &~ ALIGN;
	NEWMAT(taxa->base, alignTotal, alignVunits);
	NEWMAT(taxa->boot, alignTotal, alignVunits);
	newmem(taxa->rdgs, alignTotal);

	// Allocate space for taxa, visiteds
	newmem(taxa->taxa, taxa->nTotal);
	newmem(taxa->visit, taxa->nTotal);
	txUnvisit(taxa);

	// Allocate Scan buffer for variants, add 1 for EOS
	newmem(vbuf, taxa->nVunits+1);

	// Scan bulk of file
	for (nt = 0; nt < taxa->nTotal; nt++) {
		Taxon *tx = &taxa->taxa[nt];
		char *vb;

		if (nt < taxa->nExtant) {
			if (fscanf(fp, "%s %s", name, vbuf) != 2)
				return 0;
			vb = vbuf;
		} else {
			sprintf(name, "[%d]", nt - taxa->nExtant + 1);
			vb = (char *) 0;
		}

		tx->name = strdup(name);
		assert( tx->name );
		taxa->rdgs[nt] = txInVars(taxa, nt, vb);
		assert( taxa->rdgs[nt] == taxa->base[nt] );

		tx->nSings = 0;

		newmem(tx->isSing, taxa->nVunits);
		ZERO(tx->isSing, taxa->nVunits);
	}
	free(vbuf);

	// Allocate and initialize remaining support vectors
	newmem(taxa->type, taxa->nVunits);
	newmem(taxa->majority, taxa->nVunits);

	for (vv = 0; vv < taxa->nVunits; vv++)
		taxa->type[vv] = txVarType(taxa, vv);	// Also sets majority

	newmem(taxa->permute, taxa->nVunits);
	taxa->perturbed = NO;

	return taxa;
}

void
txFree(Taxa **taxa)
{
	Cursor v;
	Taxa *tx = *taxa;

	for (v = 0; v < tx->nTotal; v++) {
		free(tx->taxa[v].name);
	}
	free(tx->majority);
	free(tx->type);
	free(tx->permute);
	FREMAT(tx->base);
	FREMAT(tx->boot);
	free(tx->taxa);
	free(tx);
	*taxa = (Taxa *) 0;
	return;
}

Cursor
txFind(const Taxa *tx, const char *name)
{
	Cursor tt;

	if (!name)
		return -1;

	for (tt = 0; tt < tx->nTotal; tt++) {
		if (strcmp(txName(tx,tt), name) == 0)
			return tt;
	}
	return -1;
}

void
	txPermuteTaxon(Taxa *taxa, Cursor tt)
{
	Cursor vv;
	vunit *base = taxa->base[tt];
	vunit *boot = taxa->boot[tt];
	vunit *rdgs = taxa->rdgs[tt];

	if (!taxa->perturbed) {
		assert( rdgs == base );
		return;
	}
	assert( rdgs == boot );
	for (vv = 0; vv < taxa->nVunits; vv++)
		boot[vv] = base[taxa->permute[vv]];
}

Length
	txApographic(Taxa *taxa, Cursor node, Cursor outgroup)
{
	Length cost = 0;
	Cursor vv;

	if (outgroup == -1) {
		for (vv = 0; vv < taxa->nVunits; vv++)
			cost += txRdgs(taxa,node)[vv] != 1;
		return cost;
	}

	// Houston, we have an outgroup.
	for (vv = 0; vv < taxa->nVunits; vv++) {
		vunit av = txRdgs(taxa,outgroup)[vv];
		cost += txRdgs(taxa,node)[vv] != av;
	}
	return cost;
}

vunit
	txVprint(Taxa *taxa, Cursor node, Cursor vv)
{
	vunit tv = txRdgs(taxa,node)[vv];

 	return STATES[tv];
}

// Perturbation for Bootstrap and Ratchet

void
	txPerturb(Taxa *tx, CodeID *codes)
{
	Cursor tt, vv;

	// Set up permutation vector
	for (vv = 0; vv < tx->nVunits; vv++) {
		Cursor vb = (int) ((double) tx->nVunits * random() / (RAND_MAX + 1.0));
		tx->permute[vv] = vb;
	}
	tx->perturbed = YES;
	
	// Permute and point rdgs now to boot[]
	for (tt = 0; tt < tx->nTotal; tt++) {
		tx->rdgs[tt] = tx->boot[tt];
		txPermuteTaxon(tx, tt);
		codes[tt] = TaxonCode++;
	}
}

void
	txRestore(Taxa *tx, CodeID *codes)
{
	Cursor tt;

	for (tt = 0; tt < tx->nTotal; tt++) {
		tx->rdgs[tt] = tx->base[tt];
		codes[tt] = TaxonCode++;
	}
	tx->perturbed = NO;
}

