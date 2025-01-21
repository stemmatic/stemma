/*
	log.c - log results
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "netmsg.h"
#include "cache.h"
#include "boot.h"
#include "log.h"

struct Log {
	int argc;
	char **argv;

	char *inbase;
	char *outbase;
};

struct NetStats {
	struct {
		Cursor fr, to;
		double pct;
		int biggest;
	} *retLinks;
	int nEnd;
	Cursor vu;
};
typedef struct NetStats NetStats;

static void logStemma(FILE *log, Net *nt);
static void logTree(FILE *log, Net *nt, NetStats *ns);
static void logLinks(FILE *log, Net *nt);
static void logChars(FILE *log, Net *nt, Cursor node, Cursor vv);
static void logMixes(FILE *log, Net *nt, NetStats *ns);

static void logAnnote(Log *lg, Net *nt, Cursor node);
static void logApographies(Log *lg, Net *nt);
static void logTaxonHeader(FILE *fpout, Net *nt, Cursor t);
static void logTaxonApographies(FILE *fpout, Net *nt, Cursor t, FILE *fpin);
static Cursor logMRCA(Net *nt, Cursor mixedT, int *pnCAs);
static int logBPviols(Net *nt, Cursor t, Cursor mrca);

/////////////////////////////////////////////////////////////////


struct LogFiles {
	enum Logs id;
	char *ext;
	char *desc;
} LogFiles[] = {
	{ lgTX,   ".tx",   "<taxon-file> Taxon-variant matrix", },
	{ lgNO,   ".no",   "<con-file>   Stratigraphic constraints", },
	{ lgVR,   ".vr",   "<var-file>   Listing of variants", },
	{ lgSTEM, ".STEM", "             Stemma, nodes, links, and mixtures", },
	{ lgAPOS, ".APOS", "             Apographies at each node", },
	{ lgNOTE, ".NOTE", "             Annotated apparatus", },
	{ lgBOOT, ".BOOT", "BOOT=<nreps> Bootstrap percentages", },
	{ lgSTAT, ".STAT", "             Overall similarity statistics", },
	{ lgMIXS, ".MIXS", "             Mixture results", },
	{ lgVARS, ".VARS", "             Collation variants for all nodes", },
};

Log *
	logInit(int argc, char *argv[], int maxArg, char *usage)
{
	static Log lg[1];			// Not reentrant, but YNGNI
	enum Logs ll;

	// Sanity check my enum, struct array
	for (ll = 0; ll < lgMAX; ll++) {
		assert( ll == LogFiles[ll].id );
	}

	if (argc < 2 || argc > maxArg) {
		fprintf(stderr, "Usage: %s infiles %s >diags\n", argv[0], usage);

		fprintf(stderr, "Input and OUTPUT used files:\n");
		for (ll = 0; ll < lgMAX; ll++) {
			fprintf(stderr, "%.6s %s\n",
				LogFiles[ll].ext, LogFiles[ll].desc);
		}
		
		exit(-1);
		return (Log *) 0;
	}

	lg->argc = argc;
	lg->argv = argv;
	lg->inbase = argv[1];
	lg->outbase = argv[2];
	if (!lg->outbase)
		lg->outbase = lg->inbase;

	return lg;
}

FILE *
	logFile(Log *lg, enum Logs logExt)
{
	char filename[80];
	char *fn = filename;
	char *base, *mode;
	char *ext = LogFiles[logExt].ext;
	int ch;
	FILE *fp;

	// Set up output mode
	if (logExt >= lgSTEM) {
		base = lg->outbase;
		mode = "w";
	} else {
		base = lg->inbase;
		mode = "r";
	}

	// Just return extension if there's no base
	if (!base)
		strcpy(filename, ext+1);
	else {
		while ((ch = *fn++ = *base++) && ch != '.')
			;
		--fn;
		while ((*fn++ = *ext++))
			;
	}

	fp = fopen(filename, mode);
	if (!fp)
		perror(filename);
	return fp;
}

void
	logAnalysis(FILE *log, Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor vv;
	enum VarType n;
	int count[nVarType];
	char *txt[nVarType] = { "C", "S", "I", };

	fprintf(log, "Witnesses: %d; ", tx->nExtant);

	for (n = 0; n < nVarType; n++)
		count[n] = 0;
	for (vv = 0; vv < tx->nVunits; vv++)
		count[tx->type[vv]]++;

	fprintf(log, "Variants:");
	for (n = 0; n < nVarType; n++) {
		if (txt[n] > 0)
			fprintf(log, " %s=%d", txt[n], count[n]);
	}
	fprintf(log, " Total=%d; ", tx->nVunits);

	fprintf(log, "R= "CST_F" ", RetCost);
	if (nt->outgroup != -1)
		fprintf(log, "Outgroup: %s\n", txName(tx,nt->outgroup));

	if (count[Constant] > 0) {
		fprintf(log, "Constant variation units:");
		for (vv = 0; vv < tx->nVunits; vv++) {
			if (tx->type[vv] == Constant)
				fprintf(log, " %d", vv);
		}
		fprintf(log, "\n");
	}

	{
		Cache *soln;
		Length cost;

		soln = cacheNew(nt);
		cost = ntCost(nt);
		ntPropagate(nt);
		cacheSave(nt, cost, soln);
		
		fprintf(log, CST_F"-"RTC_F" u%.3f ",
			cacheCost(soln), cacheRootCost(soln),
			(double) cacheCost(soln)/cacheNodes(soln));
		fprintf(log, "x%3d\n", cacheMixedNodes(soln));
		cacheFree(nt, soln);
	}
}

static void
	logTaxonHeader(FILE *fpout, Net *nt, Cursor t)
{
	Taxa *tx = nt->taxa;
	Cursor p;

	fprintf(fpout, "Node %s (cost=" CST_F, txName(tx,t), nt->cumes[t]);
	if (nt->nParents[t] > 1)
		fprintf(fpout, ", %ld mix cost", RetCost * (nt->nParents[t]-1));
	if (txFrag(tx,t))
		fprintf(fpout, ", fragmentary");
	fprintf(fpout, "):\n");

	if (txCorrecting(tx,t) != TXNOT)
		fprintf(fpout, "Correcting %s; ", txName(tx,txCorrecting(tx,t)));

	fprintf(fpout, "Extant descendents of %s:", txName(tx,t));
	for (p = 0; p < tx->nExtant; p++) {
		if (nt->inuse[t] && nt->descendents[t][p])
			fprintf(fpout, " %s", txName(tx,p));
	}
	fprintf(fpout, "\n\n");
}


static void
	logTaxonApographies(FILE *fpout, Net *nt, Cursor t, FILE *fpin)
{
	Taxa *tx = nt->taxa;
	Cursor vv;
	char verse[256];
	char line[256];
	char buffer[256];

	int newVerse, newLine;
	Cursor p;

	logTaxonHeader(fpout, nt, t);

	rewind(fpin);
	ntCost(nt);		// Get the states
	newVerse = newLine = NO;
	for (vv = 0; vv < tx->nVunits; vv++) {
		vunit tv = txVprint(tx,t,vv);
		Cursor up, dn;
		vunit pv;
		unsigned wgt = 1;

		if (nt->vcs[t][vv] == 0 && (nt->nParents[t] > 0 || tv == '0'))
			continue;

		// Advance in vr file to variation unit
		while (fgets(buffer, sizeof buffer, fpin)) {
			Cursor v;
			char *end;
			v = (Cursor) strtol(buffer, &end, 10);

			// Weighted vunits only have the last number listed in the
			// .vr file, so advance vv (and tv) to what we've read.
			if (v > vv) {
				wgt = v-vv+1;
				vv = v;
				tv = txVprint(tx,t,vv);
			}

			if (buffer != end && v == vv)
				break;
			else if (*buffer == '@') {
				strcpy(verse, buffer);
				newVerse = YES;
			} else if (*buffer == '>') {
				strcpy(line, buffer);
				newLine = YES;
			}
		}
		
		if (newVerse) {
			fputs(verse, fpout);
			newVerse = NO;
		}
		if (newLine) {
			fputs(line, fpout);
			newLine = NO;
		}
		fprintf(fpout, "%c%c%s", tv,
			(tx->type[vv] == Singular || txIsSingular(tx,t)[vv])
				? '*' : ' ',
			buffer);

		{
			fprintf(fpout, "       ");
			for (up = nt->Ups[t]; up != -1; up = nt->br[up].nxtUp) {
				p = nt->br[up].fr;
				pv = txVprint(tx,p,vv);
				fprintf(fpout, " %s(%c)", txName(tx,p), pv);
			}
			fprintf(fpout, " -> %s(%c)", txName(tx,t), tv);
			if (nt->nChildren[t] > 0)
				fprintf(fpout, " ->");
			for (dn = nt->Dns[t]; dn != -1; dn = nt->br[dn].nxtDn) {
				p = nt->br[dn].to;
				pv = txVprint(tx,p,vv);
				fprintf(fpout, " %s(%c)", txName(tx,p), pv);
			}
			fprintf(fpout, "\n\n");
		}
		logChars(fpout, nt, t, vv);

		if (wgt > 1)
			fprintf(fpout, "Weight: %u\n\n", wgt);
	}

	fprintf(fpout, "\n");
}

static void
	logApographies(Log *lg, Net *nt)
{
	FILE *fpin, *fpout;
	Cursor t;

	fpin = logFile(lg, lgVR);
	if (!fpin) {
		perror("logApographies note-infile");
		return;
	}

	fpout = logFile(lg, lgAPOS);
	if (!fpout) {
		perror("logApographies note-outfile");
		fclose(fpin);
		return;
	}

	for (t = 0; t < nt->maxTax; t++) {
		if (!nt->inuse[t])
			continue;
		logTaxonApographies(fpout, nt, t, fpin);
	}
	fprintf(fpout, "\n");

	fclose(fpout);
	fclose(fpin);
}

#define MMAX (5) // Cheap limit

static void
	logMixParentage(FILE *log, Net *nt, Cursor t, NetStats *ns, Length *totStretch)
{
	Taxa *tx = nt->taxa;
	Cursor p1, p2, vv;
	int nn;
	double mrdgs[5];		// Mixed readings per parentage
	Cursor mpars[5];		// Identities of the parentage
	int mps = 0;			// Number of mixed parents
	double nMixRdgs = 0.0;	// Number of mixed readings
	int nMiss = 0;			// Number of missing parents states
	char *badMix = "";

#if UNMIX
		for (Cursor *un = nt->unMixed; un < nt->endMixed; un++) {
			if (!nt->descendents[t][*un])
				continue;
			badMix = "!!!";
			break;
		}
#else
	if (nt->banMixed[t])
		badMix = "!";
#endif
	fprintf(log, "%sMixed Node %s (cost "CST_F" > RetCost"CST_F"):\n",
		badMix, txName(tx,t), nt->cumes[t], RetCost*(nt->nParents[t]-1));

	// Print individual links
	for (p1 = 0; p1 < nt->maxTax; p1++) {
		if (!nt->inuse[p1])
			continue;
		if (!nt->connection[p1][t])
			continue;
		nn = 0;
		nMiss = 0;
		for (vv = 0; vv < tx->nVunits; vv++) {
			if (txRdgs(tx,p1)[vv] == MISSING && txRdgs(tx,t)[vv] != MISSING)
				nMiss++;
			if (txRdgs(tx,t)[vv] != txRdgs(tx,p1)[vv])
				continue;
			for (p2 = 0; p2 < nt->maxTax; p2++) {
				if (!nt->inuse[p2])
					continue;
				if (!nt->connection[p2][t])
					continue;
				if (p2 == p1)
					continue;
				if (txRdgs(tx,t)[vv] == txRdgs(tx,p2)[vv])
					break;
			}
			if (p2 < nt->maxTax)
				continue;
			nn++;
		}
		if (mps == MMAX)
			continue;
		nMixRdgs += nn;
		mrdgs[mps] = nn;
		mpars[mps] = p1;
		mps++;
	}

	// Calc & print the mixture percentages
	for (nn = 0; nn < mps; nn++) {
		double linkPct = 100.0*mrdgs[nn]/nMixRdgs;
		int pp = mpars[nn];
		fprintf(log, "\t"CST_F" %4s:",
			txPole(tx,pp), txName(tx,pp));
		fprintf(log, " ---> "CST_F"",
			(Length) mrdgs[nn]);
		fprintf(log, " %.0f%%\n",
			linkPct);

		// Collect net stats for identifying the biggest share
		if (ns) {
			Cursor rr;
			int amBiggest = YES;
			for (rr = 0; rr < ns->nEnd; rr++) {
				if (ns->retLinks[rr].to != t)
					continue;
				if (ns->retLinks[rr].pct < linkPct)
					ns->retLinks[rr].biggest = NO;
				else
					amBiggest = NO;
			}
			ns->retLinks[ns->nEnd].fr = pp;
			ns->retLinks[ns->nEnd].to = t;
			ns->retLinks[ns->nEnd].pct = linkPct;
			ns->retLinks[ns->nEnd].biggest = amBiggest;
			ns->nEnd++;
		}
	}

	// Common ancestors & bi-polarity violations
	block {
		int nCAs;
		int mrca = logMRCA(nt, t, &nCAs);
		int nBPviols = logBPviols(nt, t, mrca);
		fprintf(log, "Common ancestors: %d; MRCA: %s @ "CST_F"~"CST_F"; bi-polarity violations: %d\n",
			nCAs, txName(tx,mrca), txApographic(tx,mrca,nt->outgroup), txPole(tx,mrca), nBPviols);
	}

	// Calc & print pole stretch
	block {
		Length hi = 0L, lo = ~0L;
		for (nn = 0; nn < mps; nn++) {
			Length pp = txPole(tx, mpars[nn]);
			if (pp > hi) hi = pp;
			if (pp < lo) lo = pp;
		}
		char *didTC = (txPole(tx,t) < lo) ? "!TC!" : "";
		fprintf(log, "Pole Stretch: hi:"CST_F" - lo:"CST_F" = "CST_F" @ "CST_F" %s\n",
			hi, lo, hi - lo, txPole(tx, t), didTC);
		*totStretch += hi - lo;
	}
 
	fprintf(log, "\n");
}

static Cursor
	logMRCA(Net *nt, Cursor mixedT, int *pnCAs)
{
	Cursor anc, mrca;
	int nCAs = 1;

	mrca = ntRoot(nt);
	nCAs = 1;
	for (anc = 0; anc < nt->maxTax; anc++) {
		Flag isCommonAncestor = YES;
		Cursor up;
		if (!nt->inuse[anc])
			continue;
		
		// Check if every parent of mixedT is a descendent of anc
		for (up = nt->Ups[mixedT]; up != -1; up = nt->br[up].nxtUp) {
			Cursor p = nt->br[up].fr;
			if (!nt->descendents[anc][p])
				isCommonAncestor = NO;
		}

		// Check if this common ancestor is more recent than/a descendent of
		// the one we've got.
		if (isCommonAncestor) {
			if (nt->descendents[mrca][anc])
				mrca = anc;
			else if (!nt->descendents[anc][mrca])
				nCAs++;
		}
	}
	if (pnCAs)
		*pnCAs = nCAs;
	return mrca;
}

static int
	logBPviols(Net *nt, Cursor t, Cursor mrca)
{
	Taxa *tx = nt->taxa;
	int nViols = 0;
	Cursor vv;

	// Just print out mix states
	for (vv = 0; vv < tx->nVunits; vv++) {
		Cursor up, p;
		Flag hasCommonReading = YES;
		Flag hasError = NO;

		for (up = nt->Ups[t]; up != -1; up = nt->br[up].nxtUp) {
			p = nt->br[up].fr;
			if (txVprint(tx,p,vv) != txVprint(tx,mrca,vv)
			&&  txVprint(tx,p,vv) != STATES[MISSING])
				hasError = YES;
			if (nt->br[up].nxtUp != -1) {
				Cursor nxtUp = nt->br[up].nxtUp;
				Cursor nxtP = nt->br[nxtUp].fr;
				if (txVprint(tx,p,vv) != txVprint(tx,nxtP,vv))
					hasCommonReading = NO;
			}
		}
		if (hasCommonReading && hasError)
			nViols += 1;
	}
	return nViols;
}

static void
	logMixtures(Log *lg, Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor t, vv;
	Length totStretch = 0;	// Total stretch distances

	FILE *fpin, *fpout;

	fpin = logFile(lg, lgVR);
	if (!fpin) {
		perror("logMixtures note-infile");
		return;
	}

	fpout = logFile(lg, lgMIXS);
	if (!fpout) {
		perror("logMixtures note-outfile");
		fclose(fpin);
		return;
	}

	for (t = 0; t < nt->maxTax; t++) {
		Cursor up, p;

		if (!nt->inuse[t])
			continue;
		if (nt->nParents[t] < 2)
			continue;

		fprintf(fpout, "Mixed ");
		logTaxonHeader(fpout, nt, t);
		logMixParentage(fpout, nt, t, (NetStats *) 0, &totStretch);

		ntCost(nt);		// Get the states

#if 1
		// Print Interleaved mixture analysis
		block {
			Cursor mrca;
			int nCAs = 0;

			mrca = logMRCA(nt, t, &nCAs);

			// Just print out mix states
			for (vv = 0; vv < tx->nVunits; vv++) {
				Flag hasCommonReading = YES;
				Flag hasError = NO;
				Flag hasApography = YES;
				char buff[512];
				char *bp = buff;
				Cursor affinity = mrca;

				bp += sprintf(bp, "%6u ", vv);
				bp += sprintf(bp, " %c", txVprint(tx,t,vv));
				bp += sprintf(bp, " <");
				for (up = nt->Ups[t]; up != -1; up = nt->br[up].nxtUp) {
					p = nt->br[up].fr;
					bp += sprintf(bp, " %c", txVprint(tx,p,vv));
					if (txVprint(tx,p,vv) != txVprint(tx,mrca,vv) && txVprint(tx,p,vv) != STATES[MISSING]) {
						hasError = YES;
						if (txVprint(tx,p,vv) == txVprint(tx,t,vv))
							affinity = p;
					}
					if (nt->br[up].nxtUp != -1) {
						Cursor nxtUp = nt->br[up].nxtUp;
						Cursor nxtP = nt->br[nxtUp].fr;
						if (txVprint(tx,p,vv) != txVprint(tx,nxtP,vv))
							hasCommonReading = NO;
					}
					if (txVprint(tx,p,vv) == txVprint(tx,t,vv))
						hasApography = NO;
				}
				if (hasApography)
					affinity = t;
				bp += sprintf(bp, " > ");
				bp += sprintf(bp, " %c", txVprint(tx,mrca,vv));
				if (hasError || hasApography)
					fprintf(fpout, "%s %6s%s\n", buff, txName(tx, affinity), (hasCommonReading && hasError) ? " !!" : "");
			}
			fprintf(fpout, "\n\n");
		}
	
#endif
#if 0
		// Go through each parent that contributed the reading.
		for (up = nt->Ups[t]; up != -1; up = nt->br[up].nxtUp) {
			int newVerse, newLine;
			char verse[256];
			char line[256];
			char buffer[256];

			p = nt->br[up].fr;

			fprintf(fpout, "Contributions from parent %s:\n\n",
				txName(tx, p));
			rewind(fpin);
			newVerse = newLine = NO;
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit tv = txVprint(tx,t,vv);
				vunit pv = txVprint(tx,p,vv);
				Cursor u2, dn;
				int mixed_here = NO;

				if (nt->vcs[t][vv] > 0)
					continue;
				if (tv == STATES[MISSING])
					continue;
				if (tv != pv)
					continue;

				// Ignore variant units where all parents agree.
				for (u2 = nt->Ups[t]; u2 != -1; u2 = nt->br[u2].nxtUp) {
					Cursor p2 = nt->br[u2].fr;
					char pv2 = txVprint(tx,p2,vv);
					if (p2 == p)
						continue;
					if (pv2 != pv)
						mixed_here = YES;
				}
				if (!mixed_here)
					continue;
	
				// Advance in vr file to variation unit
				while (fgets(buffer, sizeof buffer, fpin)) {
					Cursor v;
					char *end;
					v = (Cursor) strtol(buffer, &end, 10);

					// Weighted vunits only have the last number listed in the
					// .vr file, so advance vv (and tv) to what we've read.
					if (v > vv) {
						vv = v;
						tv = txVprint(tx,t,vv);
					}

					if (buffer != end && v == vv)
						break;
					else if (*buffer == '@') {
						strcpy(verse, buffer);
						newVerse = YES;
					} else if (*buffer == '>') {
						strcpy(line, buffer);
						newLine = YES;
					}
				}
				
				if (newVerse) {
					fputs(verse, fpout);
					newVerse = NO;
				}
				if (newLine) {
					fputs(line, fpout);
					newLine = NO;
				}
				fprintf(fpout, "%c%c%s", tv,
					(tx->type[vv] == Singular || txIsSingular(tx,t)[vv])
						? '*' : ' ',
					buffer);

				{
					fprintf(fpout, "       ");
					for (u2 = nt->Ups[t]; u2 != -1; u2 = nt->br[u2].nxtUp) {
						Cursor p2 = nt->br[u2].fr;
						vunit pv2 = txVprint(tx,p2,vv);
						fprintf(fpout, " %s(%c)", txName(tx,p2), pv2);
					}
					fprintf(fpout, " -> %s(%c)", txName(tx,t), tv);
					if (nt->nChildren[t] > 0)
						fprintf(fpout, " ->");
					for (dn = nt->Dns[t]; dn != -1; dn = nt->br[dn].nxtDn) {
						Cursor p2 = nt->br[dn].to;
						vunit pv2 = txVprint(tx,p2,vv);
						fprintf(fpout, " %s(%c)", txName(tx,p2), pv2);
					}
					fprintf(fpout, "\n\n");
				}
			}
		}
		fprintf(fpout, "\n");
#endif
	}
	fprintf(fpout, "Total stretch: "CST_F"\n", totStretch);
	fprintf(fpout, "\n");

	fclose(fpout);
	fclose(fpin);
}

static void
	logStemma(FILE *log, Net *nt)
{
	Link *work=0, *link=0, *end=0;

	end = ntPreorderLinks(nt, &work);

	// Loop over each character
	for (link = work; link < end; link++) {
		fprintf(log, ""NAM_F" "LNK_F" "NAM_F"; ",
			txName(nt->taxa,link->from),
			C_NLINK(nt,link->to),
			txName(nt->taxa,link->to));
	}
	fprintf(log, "\n");
	free(work);
}

static double
	nsParentage(NetStats *ns, Cursor fr, Cursor to)
{
	Cursor t;
	for (t = 0; t < ns->nEnd; t++) {
		if (ns->retLinks[t].fr == fr && ns->retLinks[t].to == to)
			return ns->retLinks[t].pct;
	}
	return 100.0;
}

static int
	nsBiggest(NetStats *ns, Cursor fr, Cursor to)
{
	Cursor t;
	for (t = 0; t < ns->nEnd; t++) {
		if (ns->retLinks[t].fr == fr && ns->retLinks[t].to == to)
			return ns->retLinks[t].biggest;
	}
	return YES;
}

static void
	nsPrintState(FILE *fp, NetStats *ns, Taxa *tx, Cursor node)
{
	if (!ns || ns->vu == ERR)
		return;
	fprintf(fp, "{%c}", txVprint(tx,node,ns->vu));
}

static void
	logInorder(Net *nt, FILE *fp, Cursor node, NetStats *ns, double pct, int big)
{
	Taxa *tx = nt->taxa;
	Cursor dn;
	char prefix[16];

	if (nt->nParents[node] > 1) {
		sprintf(prefix, "%.0f%% ", pct);
	} else
		prefix[0] = EOS;
	
	if (big == NO) {
		fprintf(fp, "'%s"NAM_F"'", prefix, txName(tx,node));
		nsPrintState(fp, ns, tx, node);
		fprintf(fp, ":%lu", (unsigned long) nt->cumes[node]);
		return;
	}

	if (nt->nChildren[node] == 0) {
		fprintf(fp, "'%s"NAM_F"'", prefix, txName(tx,node));
	} else {
		fprintf(fp, "(");
		if (nt->nParents[node] > 1 || node < tx->nExtant) {
			// Remove this conditional to label every internal node.
			fprintf(fp, "'%s"NAM_F"':0, ", prefix, txName(tx,node));
		}
		for (dn = nt->Dns[node]; dn != -1; dn = nt->br[dn].nxtDn) {
			Cursor kid = nt->br[dn].to;
			pct = (nt->nParents[kid] > 1) ? nsParentage(ns, node, kid) : 100.0;
			big = (nt->nParents[kid] > 1) ? nsBiggest(ns, node, kid) : YES;
			logInorder(nt, fp, kid, ns, pct, big);
			if (nt->br[dn].nxtDn != -1)
				fprintf(fp, ", ");
		}
		fprintf(fp, ")");
	}
	nsPrintState(fp, ns, tx, node);
	fprintf(fp, ":%lu", (unsigned long) nt->cumes[node]);
	return;
}

static void
	logTree(FILE *log, Net *nt, NetStats *ns)
{
	Cursor tt;
	
	// Print a tree for every root (parentless) node.
	for (tt = 0; tt < nt->maxTax; tt++) {
		if (nt->nParents[tt] == 0 && nt->inuse[tt]) {
			logInorder(nt, log, tt, ns, 100.0, YES);
			fprintf(log, ";\n");
		}
	}
	fprintf(log, "\n");
}


static void
	logLinks(FILE *log, Net *nt)
{
	Taxa *tx = nt->taxa;
	Link *work=0, *link=0, *end=0;

	end = ntPreorderLinks(nt, &work);
	if (work == end) {
		fprintf(log, "No links!\n");
		free(work);
		return;
	}

	fprintf(log, "                               Cost  Dist~Normd Delta ?Miss\n");
	fprintf(log, "Root: "NAM_F" %*s %3ld "CST_F"~"CST_F"\n",
		txName(tx,work->from),
		24 - (int) strlen(txName(tx,work->from)), "|  ",
		nt->cumes[work->from],
		txApographic(tx,work->from,nt->outgroup),
		txPole(tx,work->from));

	// Loop over each character
	for (link = work; link < end; link++) {
		Cursor to = link->to;
		Cursor fr = link->from;
		fprintf(log, "Link: "NAM_F" "LNK_F" "NAM_F" ",
			txName(tx,fr), C_NLINK(nt,to), txName(tx,to));
		fprintf(log, "%*s", 20 - (int) strlen(txName(tx,fr)) - (int) strlen(txName(tx,to)), "|  ");
		fprintf(log, " %3ld", nt->cumes[to]);
		fprintf(log, " "CST_F"", txApographic(tx,to,nt->outgroup)); 
		fprintf(log, "~"CST_F"", txPole(tx,to));
		fprintf(log, " %5ld", txApographic(tx,to,nt->outgroup)
								- txApographic(tx,fr,nt->outgroup));
#if 1 /* was: DO_MP2 */
		{
			int vv, nMiss = 0;
			for (vv = 0; vv < tx->nVunits; vv++) {
				if (txRdgs(tx,fr)[vv] == MISSING
				&& txRdgs(tx,to)[vv] != MISSING)
					nMiss++;
			}
			if (nMiss > 0) fprintf(log, " ?%d", nMiss);
		}
#endif
		fprintf(log, "\n");
	}
	fprintf(log, "\n");
	free(work);
}

void
	logStats(FILE *log, Net *nt)
{
	Taxa *tx = nt->taxa;
	Cursor t1, t2;
	Cursor vv;
	Length diff, base;
	double agr;
	unsigned nSings = 0;

	if (!log)
		log = stderr;

	fprintf(log, "Similarities over all readings:\n");
	for (t1 = 0; t1 < nt->maxTax; t1++) {
		if (!nt->inuse[t1])
			continue;
		fprintf(log, "All 999 Node %5s:{\n", txName(tx,t1));
		for (t2 = 0; t2 < tx->nExtant; t2++) {
			if (!nt->inuse[t2])
				continue;
			fprintf(log, "%5s, ", txName(tx,t2));
			diff = base = 0;
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit v1 = txVprint(tx,t1,vv);
				vunit v2 = txVprint(tx,t2,vv);
				if (v1 == '?' || v2 == '?')
					continue;
				base++;
				if (v1 != v2)
					diff++;
			}
			agr = 0.0;
			if (base > 0)
				agr = 100.0*((double) base - diff)/(double) base;
			fprintf(log, "%.2f\n", agr);
		}
		fprintf(log, "}\n");
	}
	
	fprintf(log, "\nNumber of Singular Readings:\n");
	for (t1 = 0; t1 < nt->maxTax; t1++) {
		if (!nt->inuse[t1])
			continue;
		fprintf(log, "%5s, %d\n", txName(tx,t1), txNSings(tx,t1));
		nSings += txNSings(tx,t1);
	}
	fprintf(log, "Average number of singular readings: %g (extant), %g (all)\n",
		(double) nSings/tx->nExtant, (double) nSings/nt->nTaxa);

	fprintf(log, "\nSimilarities over informative readings:\n");
	for (t1 = 0; t1 < nt->maxTax; t1++) {
		if (!nt->inuse[t1])
			continue;
		fprintf(log, "Informative 999 Node %5s:{\n", txName(tx,t1));
		for (t2 = 0; t2 < tx->nExtant; t2++) {
			if (!nt->inuse[t2])
				continue;
			fprintf(log, "%5s, ", txName(tx,t2));
			diff = base = 0;
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit v1 = txVprint(tx,t1,vv);
				vunit v2 = txVprint(tx,t2,vv);
				if (v1 == '?' || v2 == '?')
					continue;
				// Count non-singulars
				if (txIsSingular(tx,t1)[vv] || txIsSingular(tx,t2)[vv])
					continue;
				base++;
				if (v1 != v2)
					diff++;
			}
			agr = 0.0;
			if (base > 0)
				agr = 100.0*((double) base - diff)/(double) base;
			fprintf(log, "%.2f\n", agr);
		}
		fprintf(log, "}\n");
	}

	char *stats = getenv("STATS");
	if (!stats || strcmp(stats, "CBGM") != 0)
		return;

	fprintf(log, "\nPotential Ancestors/Descendents (CBGM):\n");
	for (t1 = 0; t1 < tx->nExtant; t1++) {
		if (!nt->inuse[t1])
			continue;
		fprintf(log, "Potential Anc/Desc Node %5s:\n", txName(tx,t1));
		for (t2 = 0; t2 < nt->maxTax; t2++) {
			int w1 = 0, w2 = 0;
			if (!nt->inuse[t2])
				continue;
			fprintf(log, "%5s, ", txName(tx,t2));
			diff = base = 0;
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit v1 = txVprint(tx,t1,vv);
				vunit v2 = txVprint(tx,t2,vv);
				if (v1 == '?' || v2 == '?')
					continue;
				if (v1 == '0' && v2 != '0')
					w1++;
				if (v2 == '0' && v1 != '0')
					w2++;
				base++;
				if (v1 != v2)
					diff++;
			}
			agr = 0.0;
			if (base > 0)
				agr = 100.0*((double) base - diff)/(double) base;
			fprintf(log, "%.1f, %d %c %d\n", agr, w1, (w1 == w2) ? '=' : (w1 > w2) ? '>' : '<', w2);
		}
		fprintf(log, "\n");
	}

	fprintf(log, "\nPotential Siblings (CBGM fix):\n");
	for (t1 = 0; t1 < tx->nExtant; t1++) {
		if (!nt->inuse[t1])
			continue;
		fprintf(log, "Potential Sibling Node %5s:\n", txName(tx,t1));
		for (t2 = 0; t2 < nt->maxTax; t2++) {
			int w1 = 0;
			if (!nt->inuse[t2])
				continue;
			fprintf(log, "%5s, ", txName(tx,t2));
			diff = base = 0;
			for (vv = 0; vv < tx->nVunits; vv++) {
				vunit v1 = txVprint(tx,t1,vv);
				vunit v2 = txVprint(tx,t2,vv);
				if (v1 == '?' || v2 == '?')
					continue;
				if (v1 == v2 && v1 != '0')
					w1++;
				base++;
				if (v1 != v2)
					diff++;
			}
			agr = 0.0;
			if (base > 0)
				agr = 100.0*((double) base - diff)/(double) base;
			fprintf(log, "%.1f, %d\n", agr, w1);
		}
		fprintf(log, "\n");
	}
}

static void
	logSupporter(FILE *log, Net *nt, Taxa *tx, Cursor t)
{
	fprintf(log, " %s%s", (nt->nParents[t]>1) ? "%" : "", txName(tx,t));
}

static void
	logRegEx(FILE *log, Net *nt, Taxa *tx, Cursor t)
{
	fprintf(log, "|%s%s%s$", (nt->nParents[t]>1) ? "% " : "^", (t < tx->nExtant) ? "" : "\\", txName(tx,t));
}

static void
	logChars(FILE *log, Net *nt, Cursor node, Cursor vv)
{
	Taxa *tx = nt->taxa;
	Cursor t;
	Flag rdgs[256];
	int rdg, nStates;
	int m=0, s=0, g=0;	// For calculation of ci, ri.
	int nodeRdg = MISSING;
	int inGroup=0, outGroup=0;		// For calculation of apography index
	int tp = 0, fp = 0, fn = 0, tn = 0;	// For calculation of phi/MCC.

	fprintf(log, "Var %d: ", vv);

	/* Find out which rdgs are attested for this vunit */
	ZERO(rdgs, dimof(rdgs));
	for (t = 0; t < nt->maxTax; t++) {
		if (!nt->inuse[t])
			continue;
		rdg = txRdgs(tx,t)[vv];
		if (nt->nParents[t] == 0)
			fprintf(log, "%s:=%c%c ", txName(tx,t), txVprint(tx,t,vv),
				(txIsSingular(tx,t)[vv]) ? '*' : ' ');
		rdgs[rdg] = YES;
	}

	/* Print apographies */
	for (t = 0; t < nt->maxTax; t++) {
		vunit tv;
		if (!nt->inuse[t])
			continue;
		tv = txRdgs(tx,t)[vv];
		if (nt->vcs[t][vv] > 0) {
			fprintf(log, "%s=%c%c ", txName(tx,t), txVprint(tx,t,vv),
				(txIsSingular(tx,t)[vv]) ? '*' : ' ');
			s++;
		}
		if (t < tx->nExtant && tv != MISSING && tv != tx->majority[vv])
			g++;
	}
	fprintf(log, "\n");

	/* Print which nodes have which reading */
	if (node >= tx->nExtant)
		nodeRdg = txRdgs(tx,node)[vv];
	nStates = (tx->type[vv] == Informative) ? MAXSTATES : MISSING+1;
	for (rdg = 0; rdg < nStates; rdg++) {
		/* Skip unattested readings */
		if (rdgs[rdg] == NO)
			continue;
		
		if (rdg != MISSING)
			m++;

		fprintf(log, "    %c :=", STATES[rdg]);
		if (node == -1) {
			/* Print regex of witnesses (for easy coloring of phylotree.hyphy.org trees) */
			for (t = 0; t < nt->maxTax; t++) {
				if (!nt->inuse[t])
					continue;
				if (t >= tx->nExtant && nt->nParents[t] == 1)
					continue;
				if (txRdgs(tx,t)[vv] == rdg)
					logRegEx(log, nt, tx, t);
			}
		} else {
			/* Print two lines: first for descendents, then for non-descendents */
			for (t = 0; t < nt->maxTax; t++) {
				if (!nt->inuse[t])
					continue;
				if (t >= tx->nExtant && nt->nParents[t] == 1)
					continue;
				if (!nt->descendents[node][t])
					continue;
				if (txRdgs(tx,t)[vv] == rdg) {
					logSupporter(log, nt, tx, t);
					if (t < tx->nExtant && rdg == nodeRdg)
						inGroup++;
					if (t < tx->nExtant) {
						if (nodeRdg == MISSING || txRdgs(tx,t)[vv] == MISSING)
							;
						else if (nodeRdg == txRdgs(tx,t)[vv])
							tp++;
						else
							fn++;
					}
				}
			}
			fprintf(log, "\n        ");
			for (t = 0; t < nt->maxTax; t++) {
				if (!nt->inuse[t])
					continue;
				if (t >= tx->nExtant && nt->nParents[t] == 1)
					continue;
				if (nt->descendents[node][t])
					continue;
				if (txRdgs(tx,t)[vv] == rdg) {
					logSupporter(log, nt, tx, t);
					if (t < tx->nExtant && rdg == nodeRdg)
						outGroup++;
					if (t < tx->nExtant) {
						if (nodeRdg == MISSING || txRdgs(tx,t)[vv] == MISSING)
							;
						else if (nodeRdg == txRdgs(tx,t)[vv])
							fp++;
						else
							tn++;
					}
				}
			}
		}
		fprintf(log, "\n");
	}

	/* Print some summary statistics, consistency index and retentiion index. */
	if (tx->type[vv] == Informative) {
		--m;
		fprintf(log, "ci=%.2f, ri=%.2f", (double) m/(double) s, 
			((double) g - (double) s) / ((double) g - (double) m));
		if (nodeRdg != MISSING) {
			double a = inGroup;
			double b = outGroup + inGroup;
			fprintf(log, ", ai=%.3f, si=%.3f", (a/b), (a/b) * ((a-1)/(b-1)));

			fprintf(log, ", tp=%d, fp=%d, fn=%d, tn=%d", tp, fp, fn, tn);
			int sq = ((tp+fn)*(fp+tn)*(tp+fp)*(fn+tn));
			double phi = (sq > 0) ? (tp*tn - fp*fn) / sqrt(sq) : 0.0;
			fprintf(log, ", phi=%.3f", phi);
		}
		fprintf(log, "\n");
	}
	fprintf(log, "\n");
}

static void
	logAnnote(Log *lg, Net *nt, Cursor node)
{
	char buffer[256];
	FILE *fpin, *fpout;

	fpin = logFile(lg, lgVR);
	if (!fpin)
		return;

	if ((fpout = logFile(lg, lgNOTE))) {
		while (fgets(buffer, sizeof buffer, fpin)) {
			Cursor v;
			char *end;
			fputs(buffer, fpout);
			v = (Cursor) strtol(buffer, &end, 10);
			if (buffer != end)
				logChars(fpout, nt, node, v);
		}
		fclose(fpout);
	}

	fclose(fpin);
}

void
	logVariants(Log *lg, Net *nt)
{
	Taxa *tx = nt->taxa;
	char buffer[256];
	FILE *fpin, *fpout;

	fpin = logFile(lg, lgVR);
	if (!fpin)
		return;

	fpout = logFile(lg, lgVARS);
	if (!fpout)
		return (void) perror("VARS");

	while (fgets(buffer, sizeof buffer, fpin)) {
		Cursor vv, rdg;
		char *end;
		
		if (buffer[0] == '-')
			continue;
		vv = (Cursor) strtol(buffer, &end, 10);
		if (buffer == end) {
			fputs(buffer, fpout);
			continue;
		}

		fprintf(fpout, "^");
		fputs(buffer, fpout);

		rdg = MAXSTATES;
		while (--rdg >= 0) {
			Cursor t;
			int nMSS = 0;

			/* Skip the collation reading (0) */
			if (STATES[rdg] == '0')
				continue;

			/* Skip rdgs with no MSS */
			for (t = 0; t < nt->maxTax; t++) {
				if (!nt->inuse[t])
					continue;
				if (txRdgs(tx,t)[vv] == rdg)
					nMSS++;
			}
			if (nMSS == 0)
				continue;

			/* Print out the witnesses for each variant */
			fprintf(fpout, "=%c ", STATES[rdg]);
			for (t = 0; t < nt->maxTax; t++) {
				if (!nt->inuse[t])
					continue;
				if (txRdgs(tx,t)[vv] == rdg)
					fprintf(fpout, " %s", txName(tx,t));
			}
			fprintf(fpout, " \n");
		}
	}
	fclose(fpout);

	fclose(fpin);
}

void
	logUncollate(Log *lg, Net *nt, char *ms)
{
	char buffer[256];
	FILE *fpin, *fpout;
	Cursor node;

	fpin = logFile(lg, lgVR);
	if (!fpin)
		return;
	if ((node = txFind(nt->taxa, ms)) == TXNOT)
		return;

#if 0
	// Re-annotate the NOTES based on the node
	logAnnote(lg, nt, node);
#endif

	sprintf(buffer, "MS-%s.txt", ms);
	if ((fpout = fopen(buffer, "w")) != NULL) {
		fprintf(fpout, "(UN)COLLATION of %s\n\n", ms);
		while (fgets(buffer, sizeof buffer, fpin)) {
			Cursor v;
			char *end;
			if (buffer[0] == '-')
				continue;
			v = (Cursor) strtol(buffer, &end, 10);
			if (buffer != end) {
				int rdg = txVprint(nt->taxa, node, v);
				if (rdg == '0')
					continue;
				fprintf(fpout, "%c%c ", rdg,
					(txIsSingular(nt->taxa,node)[v]) ? '*' : ' ');
			}
			fputs(buffer, fpout);
		}
		fclose(fpout);
	}

	fclose(fpin);
}

void
	logResults(Log *lg, Net *nt, char *outbase)
{
	Length cost;
	FILE *log;
	char *envp;
	struct NetStats ns[1] = { { 0, }, };

	cost = ntCost(nt);
	ntPropagate(nt);

	block {
		Cursor t;
		char *vu = getenv("VU");
		newmem(ns->retLinks, nt->nLinks);
		for (t = 0; t < nt->nLinks; t++) {
			ns->retLinks[t].fr = ns->retLinks[t].to = TXNOT;
			ns->retLinks[t].pct = 0.0;
			ns->retLinks[t].biggest = NO;
		}
		ns->nEnd = 0;
		ns->vu = (vu) ? atoi(vu) : ERR;  // Which optional vu to annotate on the tree
	}

	lg->outbase = outbase;
	log = logFile(lg, lgSTEM);
	if (!log) {
		fprintf(stderr, "Cannot write results to file, using stderr...\n");
		log = stderr;
	}

	fprintf(log, "***********************\n");
	if (lg->argv) {
		int n;
		for (n = 0; n < lg->argc; n++)
			fprintf(log, "%s ", lg->argv[n]);
		fprintf(log, "\n");
	}
	fprintf(log, "Compile-time switches: DO_MAXMIX=%d DO_POLE=%d NO_TRIPS=%d UNMIX=%d NO_TC=%d IMIXD=%d IPOLE=%d ROOTFIX=%d\n",
		DO_MAXMIX, DO_POLE, NO_TRIPS, UNMIX, NO_TC, IMIXD, IPOLE, ROOTFIX);
	fprintf(log, "***********************\n\n");

	logAnalysis(log, nt);

	fprintf(log, "\n%d/%ld {%ld} ", nt->nLinks, cost, 0L);
	logStemma(log, nt);
	fprintf(log, "\n");

	// Print apographies by taxon
	logApographies(lg, nt);

	// Print apographic taxa by variant
	logAnnote(lg, nt, -1);

	logLinks(log, nt);
	logMixes(log, nt, ns);
	logMixtures(lg, nt);
	logTree(log, nt, ns);
	logTree(stdout, nt, ns);

	if ((envp = getenv("STATS")))
		logStats(logFile(lg, lgSTAT), nt);

	logVariants(lg, nt);

	fclose(log);
}

static void
	logMixes(FILE *log, Net *nt, NetStats *ns)
{
	Cursor t;
	int nMixes = 0;
	Length totStretch = 0;	// Total stretch distances

	for (t = 0; t < nt->maxTax; t++) {
		if (!nt->inuse[t])
			continue;
		if (nt->nParents[t] < 2)
			continue;

		nMixes++;
		logMixParentage(log, nt, t, ns, &totStretch);
	}
	fprintf(log, "Number of mixed nodes: %d\n", nMixes);
	fprintf(log, "Total stretch: "CST_F"", totStretch);
	fprintf(log, "\n\n");
}

void
	logBoot(Net *nt, Boot *bt, Log *lg)
{
	FILE *btfp;

	btfp = logFile(lg, lgBOOT);
	if (!btfp)
		return (void) perror("BOOT");

	bootPrint(btfp, nt, bt);

	fclose(btfp);
}

///////////////////// Various Utilities ////////////////////

#if 0
static double
	pear(double e0, double e1, double o0, double o1)
{
	double chi2 = 0.0, e = 0.0, o = 0.0;
	double etot = e0 + e1;
	double otot = o0 + o1;

	e = e0 / etot * otot;
	o = e - o0;
	chi2 += o*o/e;

	e = e1 / etot * otot;
	o = e - o1;
	chi2 += o*o/e;

	return chi2;
}

/* Convert Chi^2 to P */
static double
	ChiSqP(int df, double x2)
{
   /* found on a web page */
   double  p, t, k, a;

   p = exp( -x2 / 2 );

   if ((df & 1) > 0)
         p *= sqrt( (x2+x2) / (atan(1.0)*4) );

   k = df;
   while (k >= 2) {
         p *= x2 / k;
         k -= 2;
   }

   t = p;  a = df;
   while (t > (1.0E-7 * p)) {
         a += 2;
         t *= x2 / a;
         p += t;
   }

   return 1.0 - p;
}
#endif
