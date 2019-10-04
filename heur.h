
/*
	heur.h
*/

extern Length	heurRatchet(Net *nt, Cache *start, Cache *best, 
					int nCycles, int nRatchets);
extern Length	heurAnneal(Net *nt, Cache *start, Cache *best, 
					int temperature);
