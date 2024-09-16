/*
	net.h
*/

typedef struct Branch {
	Cursor nxtBr;		// Next branch;
	Cursor fr;			// Parent, from, node
	Cursor to;			// Child, to, node
	Cursor nxtUp;		// Next up/parent/from
	Cursor nxtDn;		// Next down/child/to
} Branch;

typedef struct Net {
	Taxa *taxa;			// Taxa (static information from .tx file)

	int nTaxa;			// Number of Active Taxa (extant + hypothetical)
	Cursor maxTax;		// One more than highest taxon
	Flag *inuse;		// Hypothetical in use?
	CodeID *codes;		// [T] Serial number of node value for caching

	Length *cumes;		// [T] Cumulative costs
	Byte **vcs;			// [T][V] Vunit costs for origination (apography)
	Byte **vcBase;		// [T][V] Base costs, for aliasing in memoized ntCost
#if DO_POLE
	Length *poles;		// [T] Pole positions
	int *banMixed;		// [T] Set to forbid a taxon from being mixed
#endif

	int nLinks;			// Number of links
	Flag **connection;	// [T][T] Connection Array, nTotal x nTotal

	Branch *br;			// [2T] Branches
	Cursor freeBr;		// Free branch
	Cursor *Ups;		// [T] Up branches (parents)
	Cursor *Dns;		// [T] Down branches (kids)

	int *nParents;		// [T] Number of parents
	int *nChildren;		// [T] Number of children

	Flag **noanc;		// [T][T] Constraint of no ancestors
	Flag **descendents;	// [T][T] Each Taxon's descendents

	Cursor outgroup;	// Outgroup taxon or -1 for '0'

#if DO_MAXMIX
	Cursor maxMix;		// Maximum number of mixed nodes
	Cursor nMixed;		// Number of mixed nodes
#endif
} Net;

typedef struct Link {
	Cursor from, to;
} Link;

extern Net		*ntNew(Taxa *tx);

extern int		ntConstraints(Net *nt, FILE *fcon);
extern int		ntConstrained(Net *nt, Cursor from, Cursor to);
extern int		ntHypConstrained(Net *nt);
extern void		ntPropagate(Net *nt);
#if DO_POLE
extern int		ntPoleCheck(Net *nt, Cursor to);
extern int		ntKidsPoleCheck(Net *nt, Cursor par);
extern void		ntPolesOK(Net *nt);
#endif

extern void		ntConnect(Net *nt, Link *link);
extern void		ntDisconnect(Net *nt, Link *link);
extern void		ntDisconnectAll(Net *nt);
extern Cursor	ntRoot(Net *nt);

extern Length	ntCost(Net *nt);
extern Length	ntKidCosts(Net *nt, Cursor node);

extern Length	ntLink(Net *nt, Length bound, Link *add);
extern Length	ntIncLink(Net *nt, Length basecost, Length bound,
					Cursor dst, Link *add);

extern Link *	ntPreorderLinks(Net *nt, Link **work);
extern void		ntRestoreLinks(Net *nt, Link *work, Link *end);
extern Cursor *	ntPostorderNodes(Net *nt, Cursor **work);
extern int		ntCyclic(Net *nt);

