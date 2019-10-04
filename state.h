
/*
	state.h
*/

extern void   stFixNode(Net *nt, Cursor node);

extern Length stAncCost(Net *nt, Cursor hyp,
				  Length *cumes, Length cost, Length bound);
extern Length stPolyCost(Net *nt, Cursor hyp, Cursor tt, Cursor t1, Cursor t2,
                  Length *cumes, Length cost, Length bound);
#ifdef ROOTFIX
extern Length stPolyRootCost(Net *nt, Cursor poly, Cursor parent, Cursor root,
                  vunit *vRoot, Length rootCost);
#endif
extern int    stLinkCost(Taxa *tx, Cursor fr, Cursor to, Flag *vc);
extern int    stRetLinkCost(Taxa *tx, Cursor fr, Cursor to, Flag *vc);
