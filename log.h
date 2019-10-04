
/*
	log.h - Log results
*/

typedef struct Log Log;

enum Logs {
	lgTAXA, lgTIME, lgVARS,
	lgSTEM, lgAPOS, lgNOTE,
	lgBOOT, lgSTAT, lgMIXS,
    lgMAX,
};

extern void logAnalysis(FILE *, Net *nt);
extern void logStats(FILE *fp, Net *nt);

extern Log *logInit(int argc, char *argv[], int maxArg, char *usage);
extern FILE *logFile(Log *lg, enum Logs logExt);
extern void logResults(Log *lg, Net *nt);
extern void logUncollate(Log *lg, Net *nt, char *ms);
extern void logBoot(Net *nt, Boot *bt, Log *lg);

extern double logVariance(FILE *log, Net *nt);
