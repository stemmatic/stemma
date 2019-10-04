
/*
	hyp.h
*/

extern Cursor   hypNew(Net *nt);	// Generate new hypothetical node
extern void     hypClear(Net *nt, Cursor hyp);    // Release hypothetical

extern void		hypFix(Net *nt);	// Fix all nodes in net

// Optimization routines
extern Length	hypLinkUp(Net *nt, Cache *curr, Cache *best);
extern Length	hypCRR(Net *nt, Cache *start, Cache *best, char *prefix);
