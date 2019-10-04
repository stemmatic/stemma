/*
	netmsg.c - Vararg messaging functions.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

#include "stemma.h"
#include "taxon.h"
#include "net.h"
#include "netmsg.h"

void
	ntVFMsg(FILE *fp, Net *nt, char *fmt, va_list ap)
{
	int c;
	static char fmtBuf[80];
	int fmtg = NO;
	char *f = fmtBuf;

	while ((c = *fmt++)) {
		switch (c) {
		case '%':
			fmtg = !fmtg;
			if (fmtg)
				f = fmtBuf;
			break;
		case 'd':
		case 'u':
		case 'c':
			if (fmtg) {
				int arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, int);
				fprintf(fp, fmtBuf, arg);
				fmtg = NO;
				continue;
			}
			break;
		case 'e':
		case 'f':
		case 'g':
			if (fmtg) {
				double arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, double);
				fprintf(fp, fmtBuf, arg);
				fmtg = NO;
				continue;
			}
			break;
		case 's':
			if (fmtg) {
				char *arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, char *);
				fprintf(fp, fmtBuf, arg);
				fmtg = NO;
				continue;
			}
			break;
		case 'N':
			if (fmtg) {
				int arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, int);
				fprintf(fp, NAM_F,
					(arg == -1) ? "(none)" : txName(nt->taxa,arg));
				fmtg = NO;
				continue;
			}
			break;
		case 'C':
			if (fmtg) {
				Length arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, Length);
				if (arg == -1)
					fprintf(fp, "MAXX");
				else
					fprintf(fp, CST_F, arg);
				fmtg = NO;
				continue;
			}
		case 'L':
			if (fmtg) {
				Link *arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, Link *);
				fprintf(fp, NAM_F LNK_F NAM_F,
					txName(nt->taxa,arg->from),
					C_NLINK(nt,arg->to),
					txName(nt->taxa,arg->to));
				fmtg = NO;
				continue;
			}
			break;
		case 'R':
			if (fmtg) {
				Length arg;
				*f++ = c;
				*f = EOS;
				arg = va_arg(ap, Length);
				if (arg == -1)
					fprintf(fp, "MX");
				else
					fprintf(fp, RTC_F, arg);
				fmtg = NO;
				continue;
			}
			break;
		default:
			break;
		}
		if (fmtg)
			*f++ = c;
		else
			putchar(c);
	}
}

void
	ntMsg(Net *nt, char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	ntVFMsg(stdout, nt, fmt, ap);
	va_end(ap);
}

void
	ntFMsg(FILE *fp, Net *nt, char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	ntVFMsg(fp, nt, fmt, ap);
	va_end(ap);
}
