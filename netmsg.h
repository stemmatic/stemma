/*
	netmsg.h
*/

// Printf-style formats for implementation

#define NAM_F "%s"
#define CST_F "%4lu"
#define RTC_F "%02lu"
#define LNK_F "%c>"
#define C_NLINK(nt,node) ((nt)->nParents[node] > 1 ? '=' : '-')

extern void     ntMsg(Net *nt, char *fmt, ...);
extern void		ntFMsg(FILE *fp, Net *nt, char *fmt, ...);
extern void		ntVFMsg(FILE *fp, Net *nt, char *fmt, va_list ap);
