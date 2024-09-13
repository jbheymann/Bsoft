/**
@file	rwIBW.h
@brief	Header file for reading and writing IGOR PRO wave files
@author Bernard Heymann
@date	Created: 20230426
@date	Modified: 20230427

	Format: Multi-dimension image file format for the IGOR PRO package
**/

#include "rwimg.h"

//#define IBWSIZE	2048

struct IBWhead {
	short	version;			// Version number for backwards compatibility.
	short	checksum;			// Checksum over this header and the wave header.
	int		wfmSize;			// The size of the WaveHeader5 data structure plus the wave data.
	int		formulaSize;		// The size of the dependency formula, including the null terminator, if any. Zero if no dependency formula.
	int		noteSize;			// The size of the note text.
	int		dataEUnitsSize;		// The size of optional extended data units.
	int		dimEUnitsSize[4];	// The size of optional extended dimension units.
	int		dimLabelsSize[4];	// The size of optional dimension labels.
	int		sIndicesSize;		// The size of string indices if this is a text wave.
	int		optionsSize1;		// Reserved. Write zero. Ignore on read.
	int		optionsSize2;		// Reserved. Write zero. Ignore on read.
} ;

struct W5head {
	int		next;				// 64	link to next wave in linked list.
	int		creationDate;		// 68	DateTime of creation.
	int		modDate;			// 72	DateTime of last modification.
	int		npnts;				// 76	Total number of points (multiply dimensions up to first zero).
	short	type;				// 80	See types (e.g. NT_FP64) above. Zero for text waves.
	short	dLock;				// 82	Reserved. Write zero. Ignore on read.
	char	whpad1[6];			// 84	Reserved. Write zero. Ignore on read.
	short	whVersion;			// 90	Write 1. Ignore on read.
	char	bname[32];			// 92	Name of wave plus trailing null.
	int		whpad2;				// 124	Reserved. Write zero. Ignore on read.
	char	dFolder[4];			// 128	Used in memory only. Write zero. Ignore on read.
	// Dimensioning info. [0] == rows, [1] == cols etc
	int		nDim[4];			// 132	Number of items in a dimension -- 0 means no data.
//	double	sfA[4];				// 148	Index value for element e of dimension d = sfA[d]*e + sfB[d].
//	double	sfB[4];
	char	sf[64];				// Dealing with non-aligned doubles
	// SI units
	char	dataUnits[4];		// 212	Natural data units go here - null if none.
	char	dimUnits[16];		// 216	Natural dimension units go here - null if none.
	// Other
	short	fsValid;			// 232	TRUE if full scale values have meaning.
	short	whpad3;				// 234	Reserved. Write zero. Ignore on read.
//	double	botFullScale;		// 236
//	double	topFullScale;		// 244	The min and max full scale value for wave.
	char	fullScale[16];		// Dealing with non-aligned doubles
	char	dataEUnits[4];		// 252	Used in memory only. Write zero. Ignore on read.
	char	dimEUnits[16];		// 256	Used in memory only. Write zero. Ignore on read.
	char	dimLabels[16];		// 272	Used in memory only. Write zero. Ignore on read.
	char	waveNoteH[4];		// 288	Used in memory only. Write zero. Ignore on read.
	char	platform;			// 292	0=unspecified, 1=Macintosh, '2=Windows', Added for Igor Pro 5.5.
	char	spare[3];			// 293
	int		whUnused[13];		// 296	Reserved. Write zero. Ignore on read.
	int		vRefNum;			// 348
	int		dirID;				// 352	Used in memory only. Write zero. Ignore on read.
	// The following stuff is considered private to Igor.
	char	priv[28];			// 356
} ;

// I/O prototypes
int 		readIBW(Bimage* p, int readdata);
int			writeIBW(Bimage* p);

