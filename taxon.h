/*
	taxon.h - Defines the Taxon Type
*/

// Valid characters in a Taxon Name.  Pretty much everything but: ,:()
#define TAXCHARS "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!@$%^&-_=+./?<>*"

// Maximum length of taxon name, mainly for I/O buffers, else dynamic
#define TAXNAME 64
#define TXNOT (-1)
#define STATES "?0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-"
#define MAXSTATES (sizeof(STATES) - 1)

// Type for variation units (i.e. characters, a term already used in C).
typedef unsigned char vunit;
#define MISSING 0

enum VarType { Constant, Singular, Informative, nVarType, };

typedef unsigned long CodeID;
extern CodeID TaxonCode;	// Used to save copy costs in caching

// Basic Information of a taxon: name and character states
typedef struct Taxon {
	char *name;
	int frag;			// Fragmentary?
	Flag *isSing;		// [nVunits] Whether the base variant is singular
	unsigned nSings;	// Number of Singular readings
	int correcting;     // Index of taxon that this one may be correcting.
} Taxon;

typedef struct Taxa {
	int nTotal;			// Number of allocated Taxa
	int nExtant;		// Number of attested, extant Taxa
	Taxon *taxa;		// [nTotal] The actual Taxa
	Flag *visit;		// [nTotal] Visited array, broken out for easy init

	int nVunits;		// Number of variation units
	enum VarType *type;	// [nVunits] Type of vunit: constant, singular, etc.
	vunit *majority;	// [nVunits] Majority character (for calc of r.i.)

	vunit **base;		// [nTotal][nVunits] Base matrix of readings
	vunit **boot;		// [nTotal][nVunits] Resampled matrix for bootstrap
	vunit **rdgs;		// [nTotal][nVunits] Points to base or boot

	int perturbed;		// Perturbed by the permutation vector?
	Cursor *permute;	// Permutation vector
} Taxa;

extern Taxa *	txScan(FILE *fp);
extern Cursor	txFind(const Taxa *taxa, const char *name);
extern void		txFree(Taxa **taxa);

extern void     txPerturb(Taxa *tx, CodeID *codes);
extern void		txRestore(Taxa *tx, CodeID *codes);
extern void		txPermuteTaxon(Taxa *taxa, Cursor node);

extern Length	txApographic(Taxa *taxa, Cursor node, Cursor outgroup);
#if DO_POLE
extern Length	txPole(Taxa *taxa, Cursor node);
#endif

extern vunit	txVprint(Taxa *taxa, Cursor node, Cursor vv);

#define txName(tx,i) (tx)->taxa[i].name
#define txFrag(tx,i) (tx)->taxa[i].frag
#define txNSings(tx,i) (tx)->taxa[i].nSings
#define txIsSingular(tx,i) (tx)->taxa[i].isSing
#define txCorrecting(tx,i) (tx)->taxa[i].correcting

#define txVisit(tx,i) (tx)->visit[i]
#define txUnvisit(tx) ZERO((tx)->visit, tx->nTotal)

#define txRdgs(tx,i) (tx)->rdgs[i]
#define txBase(tx,i) (tx)->base[i]
