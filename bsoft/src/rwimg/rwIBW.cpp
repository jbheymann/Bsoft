/**
@file	rwIBW.cpp
@brief	Functions for reading and writing image plate reader files
@author Bernard Heymann
@date	Created: 20230426
@date 	Modified: 20230427
**/

#include "rwIBW.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an IGOR PRO file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).
An image format used in the IGOR PRO package.
	File format extensions:  	.IBW
	Byte order:				flagged.
	Data type: 				various.
**/
int 	readIBW(Bimage* p, int readdata)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIBW: header sizes: " << sizeof(IBWhead) << tab << sizeof(W5head) << endl;
	
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIBW: File opened" << endl;

	IBWhead*		header = new IBWhead;
	
	fimg->read((char *)header, sizeof(IBWhead));
	if ( fimg->fail() ) return -2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIBW: Binary header read" << endl;

	W5head*			w5head = new W5head;
	
	fimg->read((char *)w5head, sizeof(W5head));
	if ( fimg->fail() ) return -3;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIBW: W5 header read" << endl;

	p->size(w5head->nDim[0], w5head->nDim[1], 1);
	p->images(w5head->nDim[2]);
	p->channels(1);

	switch ( w5head->type ) {
		case 0: p->data_type(Bit); break;
		case 2: p->data_type(Float); break;
		case 3: p->data_type(Float); p->compound_type(TComplex); p->channels(2); break;
		case 4: p->data_type(Double); break;
		case 5: p->data_type(Double); p->compound_type(TComplex); p->channels(2); break;
		case 8: p->data_type(SCharacter); break;
		case 16: p->data_type(Short); break;
		case 32: p->data_type(Integer); break;
		case 72: p->data_type(UCharacter); break;
		case 80: p->data_type(UShort); break;
		case 96: p->data_type(UInteger); break;
		default: p->data_type(UCharacter); break;
	}

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readIBW: creationDate = " << w5head->creationDate << endl;
		cout << "DEBUG readIBW: modDate = " << w5head->modDate << endl;
		cout << "DEBUG readIBW: npnts = " << w5head->npnts << endl;
//		cout << "DEBUG readIBW: sfA = " << w5head->sfA[0] << tab << w5head->sfA[1] << endl;
//		cout << "DEBUG readIBW: sfB = " << w5head->sfB[0] << tab << w5head->sfB[1] << endl;
		cout << "DEBUG readIBW: platform = " << w5head->platform << endl;
	}
	
	p->data_offset(sizeof(IBWhead) + sizeof(W5head));

	p->label(w5head->bname);
	
	double*		d = (double *)w5head->sf;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readIBW: sf = ";
		for ( int i=0; i<8; ++i ) cout << tab << d[i];
		cout << endl;
	}

	Vector3<double>		sam(d[0],d[1],1), ori(d[4],d[5],0);
	if ( strncmp(w5head->dimUnits, "m", 1) == 0 ) {
		sam *= 1e10;
		ori *= 1e10;
	}
	sam[2] = 1;
	
	p->sampling(sam);
	p->origin(ori);
	p->unit_cell_default();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readIBW: unit cell = " << p->unit_cell() << endl;
	

	if ( readdata ) {
		p->data_alloc();
		fimg->read((char *)p->data_pointer(), p->alloc_size());
		if ( fimg->fail() ) return -3;
//		if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	
	delete fimg;
	delete header;
	delete w5head;
		  
	return 0;
}

/**
@brief	Writing  an IGOR PRO file format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).
**/
int			writeIBW(Bimage* p)
{
	IBWhead*	header = new IBWhead;
	memset(header, 0, sizeof(IBWhead));
	
	header->version = 5;
	header->wfmSize = sizeof(W5head) + p->alloc_size();
	
	W5head*		w5head = new W5head;
	memset(w5head, 0, sizeof(W5head));

	w5head->npnts = p->image_size()*p->images();

	w5head->nDim[0] = p->sizeX();
	w5head->nDim[1] = p->sizeY();
	w5head->nDim[2] = p->images();
	
	double*		d = (double *)w5head->sf;
	d[0] = 1e-10*p->image->sampling()[0];
	d[1] = 1e-10*p->image->sampling()[1];
	d[2] = 1;
	d[3] = 1;
	d[4] = 0;
	d[5] = 0;
	d[6] = 0;
	d[7] = 0;
	
	w5head->dataUnits[0] = 'm';
	w5head->dimUnits[0] = 'm';
	w5head->dimUnits[4] = 'm';

	p->label().copy(w5head->bname, 31);

	p->data_offset(sizeof(IBWhead) + sizeof(W5head));

	switch ( p->data_type() ) {
		case UCharacter: w5head->type = 72; break;
		case SCharacter: w5head->type = 8; break;
		case UShort: w5head->type = 80; break;
		case Short: w5head->type = 16; break;
		case UInteger: w5head->type = 96; break;
		case Integer: w5head->type = 32; break;
		case Float: w5head->type = 2; if ( p->compound_type() == TComplex ) w5head->type = 3; break;
		case Double: w5head->type = 4; if ( p->compound_type() == TComplex ) w5head->type = 5; break;
		default: w5head->type = 72; break;
	}

	ofstream		fimg(p->file_name());
	if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, sizeof(IBWhead));
	fimg.write((char *)w5head, sizeof(W5head));
	if ( p->data_pointer() ) fimg.write((char *)p->data_pointer(), p->alloc_size());
	
	fimg.close();
	
	delete header;
	delete w5head;
		
	return 0;
}

