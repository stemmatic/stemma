
/*
	cache.h - Solution cache/checkpoint
*/

typedef struct Cache Cache;			// Opaque type

// Detailed Output
enum verbosity {
	C_VNONE,			// 0 = no output
	C_VRUN,				// 1 = run level (e.g. nTaxa stuff)
	C_VPASS,			// 2 = each pass of doing the ratchet in Anneal
	C_VITER,			// 3 = iteration level (each ratchet step)
	C_VSWAP,			// 4 = swap level: score at each step
	C_VCACHE,			// 5 = cache level: show what gets cached
	C_VLINK,			// 6 = link level: show links being considered
	C_VMAX,				// 7 = max verbosity level (default)
};

extern Cache *	cacheNew(Net *nt);
extern void		cacheFree(Net *nt, Cache *cache);

extern int		cacheBetter(Net *nt, Length cost, Length rootCost,
					int nLinks, Cache *cache);
extern int		cacheComp(Net *nt, Cache *a, Cache *b);
extern int		cacheSave(Net *nt, Length cost, Cache *cache);
extern Length	cacheRestore(Net *nt, Cache *cache);
extern void		cacheNodeRestore(Net *nt, Cursor node, Cache *cache);

extern Link *	cacheListLinks(Cache *cache, Link **work);

extern void		cacheReset(Cache *cache);
extern void		cacheMark(Cache *cache);
extern int		cacheCached(Cache *cache);

extern Length	cacheCost(Cache *cache);
extern int		cacheLinks(Cache *cache);
extern int		cacheNodes(Cache *cache);
extern Length	cacheRootCost(Cache *cache);
extern int		cacheMixedNodes(Cache *cache);

extern enum verbosity cacheVerbosity(Cache *cache);
extern void cacheSetVerbosity(Cache *cache, enum verbosity);
extern void cacheMsg(Net *nt, Cache *cache, enum verbosity, char *fmt, ...);

extern int cacheOpen(Cache *cache, char *name);
extern int cacheLock(Cache *cache);
extern int cacheUnlock(Cache *cache);
extern int cacheClose(Cache *cache);

extern int cacheWrite(Net *nt, Cache *cache);
// Return 1 if a better cached soln is read, 0 if a tie (and not read), -1 if the cached is worse.
extern int cacheRead(Net *nt, Cache *cache);
