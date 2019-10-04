/*
	stemma.h - General types, etc.
*/

// Compile-time switches for various reticulation-related issues
#define DO_MP2 0
#define ROOTFIX 0

// CLINK - Constant link (as between an original and a corrector)
#define CLINK 0
#define CLINK_FR 0
#define CLINK_TO 1

#if DO_MP2
#define MP2(a,b) (b)
#else
#define MP2(a,b) (a)
#endif

#define MMP(a,b) (a)
#define RRP(a,b) (a)

#define MAGIC (0xC0DE000A)

// Some types in all modules
typedef short Cursor;
typedef unsigned char Flag;
typedef unsigned char Byte;
typedef unsigned char Stratum;


// Tree length type
typedef unsigned long Length;
extern Length RetCost;
extern Length MaxMP2;

// Some standard CCI/ICLisms:
#define NO  0
#define YES 1
#define EOS '\0'
#define OK 0
#define ERR (-1)

#define dimof(a) ((sizeof a)/(sizeof a[0]))

#define ASET(a,c,n) memset((a), (c), (n)*sizeof *(a))
#define ACPY(a,b,n) memcpy((a), (b), (n)*sizeof *(a))
#define ZERO(a,n)   memset((a), 0x0, (n)*sizeof *(a))

// Memory management
#define newmem(p, n) if (!((p) = malloc((n) * sizeof (p)[0]))) { \
	fprintf(stderr, "No memory for %s by %lu in %s:%d (aborting).\n", \
		#p, (unsigned long) n, __FILE__, __LINE__); \
	abort(); \
} else

#define NEWMEM(p, n) if (!p && !((p) = malloc((n) * sizeof (p)[0]))) { \
	fprintf(stderr, "No memory for %s by %lu in %s:%d (aborting).\n", \
		#p, (unsigned long) n, __FILE__, __LINE__); \
	abort(); \
} else

#define NEWMAT(a, x, y) do { \
	int _ii; \
	newmem(a, (x)); \
	newmem(a[0], (x)*(y)); \
	for (_ii = 1; _ii < (x); _ii++) \
		a[_ii] = a[0] + _ii*(y); \
	} while (0)

#define FREMAT(a) (void)(free(a[0]),free(a))

// Used for moving ulong, not byte data through the bus.
#define ALIGN (sizeof (unsigned long) - 1)

#define block switch(0) default:
