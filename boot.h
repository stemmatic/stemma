
/*
	boot.h - for the bootstrap and other resampling stats
*/

// Opaque type
typedef struct Boot Boot;

extern Boot *	bootNew(Net *nt, Cache *soln);
extern void		bootStrap(Net *nt, Boot *bt, int nReps);
extern void		bootCount(Net *nt, Boot *bt);
extern void		bootPrint(FILE *, Net *nt, Boot *bt);
extern void		bootFree(Net *nt, Boot *bt);
