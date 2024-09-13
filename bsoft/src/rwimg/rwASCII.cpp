/**
@file	rwASCII.cpp
@brief	Functions for reading and writing ASCII files
@author Bernard Heymann
@date	Created: 20000318
@date	Modified: 20230526
**/

#include "rwASCII.h"
#include "Complex.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an ASCII or text image format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).

This function reads an ASCII image file with up to five dimensions
	in the order:
			c (channels), x, y, z, n (number of images)
	Default data type is Float
	Data is given as real (R) or complex (R and I or A and P),
			and may include an optional FOM (F)
	Column labels:
			Images: C X Y Z N R I F
			Structure factors: H K L A P R I F
**/
int 	readASCII(Bimage* p, int readdata)
{
    ifstream		fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long 			i, j, k, l, m(1);
	vector<long>	order(16,-1);
	vector<double>	v(16,0), min(16,0), max(16,0);
    string			s;
	Vector3<double>	ori;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: Start reading text image" << endl;
	
	// The first line contains the column labels
	getline(fimg, s);
    for ( size_t i=0; i<s.size(); i++ ) {
    	if ( isspace(s[i]) || ispunct(s[i]) ) s[i]=' ';
    	s[i] = toupper(s[i]);
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: Header line: " << s << endl;

	vector<string>	label = split(s);
	m = label.size();

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readASCII: " << m << " column labels:";
		for ( auto s: label ) cout << s << " ";
		cout << endl;
	}

	if ( verbose & VERB_DEBUG ) {
		if ( p->fourier_type() == NoTransform ) cout << "DEBUG readASCII: An image" << endl;
		else cout << "DEBUG readASCII: A transform" << endl;
	}
	
	// Associate each label with an index and a value variable
	for ( i=0; i<m; ++i ) {
		s = label[i];
		if ( order[0] < 0 && s[0] == 'X' ) order[0] = i;
		if ( order[0] < 0 && s[0] == 'H' ) {
			order[0] = i;
			p->compound_type(TComplex);
			p->fourier_type(Standard);
			p->channels(2);
		}
		if ( order[1] < 0 && s[0] == 'Y' ) order[1] = i;
		if ( order[1] < 0 && s[0] == 'K' ) order[1] = i;
		if ( order[2] < 0 && s[0] == 'Z' ) order[2] = i;
		if ( order[2] < 0 && s[0] == 'L' ) order[2] = i;
		if ( order[3] < 0 && s[0] == 'N' ) order[3] = i;
		if ( order[4] < 0 && s[0] == 'R' ) order[4] = i;
		if ( order[4] < 0 && s[0] == 'A' ) order[4] = i;
		if ( order[4] < 0 && s == "VX" ) {
			order[4] = i;
			p->compound_type(TVector3);
			p->channels(3);
		}
		if ( order[4] < 0 && s == "Red" ) {
			order[4] = i;
			p->compound_type(TRGB);
			p->channels(3);
		}
		if ( order[4] < 0 && s == "Cyn" ) {
			order[4] = i;
			p->compound_type(TCMYK);
			p->channels(4);
		}
		if ( order[5] < 0 && s[0] == 'I' ) order[5] = i;
		if ( order[5] < 0 && s[0] == 'P' ) order[5] = i;
		if ( order[5] < 0 && s == "VY" ) order[5] = i;
		if ( order[5] < 0 && s == "Grn" ) order[5] = i;
		if ( order[5] < 0 && s == "Mag" ) order[5] = i;
		if ( order[6] < 0 && s[0] == 'P' ) order[6] = i;
		if ( order[6] < 0 && s == "VZ" ) order[6] = i;
		if ( order[6] < 0 && s == "Blu" ) order[6] = i;
		if ( order[6] < 0 && s == "Yel" ) order[6] = i;
		if ( order[7] < 0 && s[0] == 'P' ) order[7] = i;
		if ( order[7] < 0 && s == "Alf" ) {
			order[7] = i;
			p->compound_type(TRGBA);
			p->channels(4);
		}
		if ( order[7] < 0 && s == "blK" ) order[7] = i;
		if ( order[8] < 0 && s[0] == 'F' ) order[8] = i;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: order: " << order[0] << " " << order[1] << " " << order[2] << " " << 
			order[3] << " " << order[4] << " " << order[5] << " " << order[6] << " " << order[7] << " " << order[8] << " (m=" << m << ")" << endl;
	
	// Pass once through the file to get the number of data points and extrema
	vector<string>	vs;
	i = 0;
 	while ( !fimg.eof() ) {
		getline(fimg, s);
		vs = split(s);
		if ( vs.size() < m ) {
			cout << "Missing data at line " << i << endl;
		} else {
			for ( j=0; j<m; ++j ) {
				v[j] = to_real(vs[j]);
				if ( min[j] > v[j] ) min[j] = v[j];
				if ( max[j] < v[j] ) max[j] = v[j];
			}
			i++;
		}
	}
	
	fimg.close();
	
	// Set some parameters for the image
	p->data_type(Float);
	p->channels(1);
	if ( order[3] > -1 ) p->images((long) (max[order[3]] + 1));
	else p->images(1);
	
	if ( p->fourier_type() == NoTransform ) {
		if ( order[0] > -1 ) {
			p->sizeX((long) (max[order[0]] - min[order[0]] + 1));
			ori[0] = min[order[0]];
		}
		if ( order[1] > -1 ) {
			p->sizeY((long) (max[order[1]] - min[order[1]] + 1));
			ori[1] = min[order[1]];
		}
		if ( order[2] > -1 ) {
			p->sizeZ((long) (max[order[2]] - min[order[2]] + 1));
			ori[2] = min[order[2]];
		}
		p->origin(ori);
	} else {
		if ( order[0] > -1 ) {
			p->sizeX((int) fabs(min[order[0]]));
			if ( p->sizeX() < fabs(max[order[0]]) ) p->sizeX((long) fabs(max[order[0]]));
			p->sizeX(p->sizeX() + 1);
		}
		if ( order[1] > -1 ) {
			p->sizeY((int) fabs(min[order[1]]));
			if ( p->sizeY() < fabs(max[order[1]]) ) p->sizeY((long) fabs(max[order[1]]));
			p->sizeY(p->sizeY() + 1);
		}
		if ( order[2] > -1 ) {
			p->sizeZ((int) fabs(min[order[2]]));
			if ( p->sizeZ() < fabs(max[order[2]]) ) p->sizeZ((long) fabs(max[order[2]]));
			p->sizeZ(p->sizeZ() + 1);
		}
	}
	
	for ( i=4; i<8; i++ ) {
		if ( order[i] > -1 ) {
			if ( p->minimum() > min[order[i]] ) p->minimum(min[order[i]]);
			if ( p->maximum() < max[order[i]] ) p->maximum(max[order[i]]);
		}
	}
	
	p->size(p->size().max(1));

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readASCII: size=" << p->size() << " channels=" << p->channels() << " images=" << p->images() << endl;

	if ( !readdata ) return 0;
	
	// Allocate memory for the data
	p->data_alloc();
	float*			fom = NULL;
	if ( order[8] > -1 ) {
		p->next = new Bimage(Float, TSimple, p->size(), p->images());
		fom = (float *) p->next->data_pointer();
	}

	// Rewind and read the data into the image
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;

	getline(fimg, s);	// Header labels

	long			c(0), x(0), y(0), z(0), n(0);
	i = 0;
	while ( !fimg.eof() ) {
		getline(fimg, s);
		vs = split(s);
		if ( vs.size() < m ) {
			cout << "Missing data at line " << i << endl;
		} else {
			for ( j=0; j<m; ++j ) v[j] = to_real(vs[j]);
			if ( order[0] > -1 ) x = (long) v[order[0]];
			if ( order[1] > -1 ) y = (long) v[order[1]];
			if ( order[2] > -1 ) z = (long) v[order[2]];
			if ( order[3] > -1 ) n = (long) v[order[3]];
			if ( x < 0 ) x += p->sizeX(); 	// Wrapping for structure factors
			if ( y < 0 ) y += p->sizeY();
			if ( z < 0 ) z += p->sizeZ();
			j = ((n*p->sizeZ() + z)*p->sizeY() + y)*p->sizeX() + x;
			for ( c=0, k=j*p->channels(), l=4; c<p->channels(); c++, k++, l++ ) p->set(k, v[order[l]]);
			if ( order[8] > -1 ) fom[j] = v[order[8]];
			i++;
			if ( verbose & VERB_DEBUG ) {
				cout << "DEBUG readASCII: " << i;
				for ( c=0, k=j*p->channels(); c<p->channels(); c++, k++ ) cout << " " << (*p)[k];
				cout << endl;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Data points read:               " << i << endl << endl;
	
	fimg.close();
	
	return 0;
}

/**
@brief	Writing an ASCII or text image format.
@param	*p			the image structure.
@return	int			error code (<0 means failure).

This function writes an ascii image file with up to five dimensions
	in the order:
			c (channels), x, y, z, n (number of images)
	Default data type is Float
	Data is given as real (R) or complex (R and I or A and P),
			and may include an optional FOM (F)
	Column labels:
			Images: C X Y Z N R I F
			Structure factors: H K L A P R I F
**/
int 	writeASCII(Bimage* p)
{
    ofstream        fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long 			i(0), nlab(128);
	long 			order[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

	string			label(nlab, ' ');
	
	// Get the data order and labels
	i = 0;
	if ( p->sizeX() > 1 ) {
		order[0] = i;
		label.replace(10*i, 1, "X");
		if ( p->fourier_type() ) label.replace(10*i, 1, "H");
		i++;
	}
	if ( p->sizeY() > 1 ) {
		order[1] = i;
		label.replace(10*i, 1, "Y");
		if ( p->fourier_type() ) label.replace(10*i, 1, "K");
		i++;
	}
	if ( p->sizeZ() > 1 ) {
		order[2] = i;
		label.replace(10*i, 1, "Z");
		if ( p->fourier_type() ) label.replace(10*i, 1, "L");
		i++;
	}
	if ( p->images() > 1 ) {
		order[3] = i;
		label.replace(10*i, 6, "Nimage");
		i++;
	}
	
	switch ( p->compound_type() ) {
		case TSimple: label.replace(10*i, 4, "Real"); i++; break;
		case TComplex: label.replace(10*i, 4, "Real"); i++;
			label.replace(10*i, 4, "Imag"); i++; break;
		case TVector2: label.replace(10*i, 2, "VX"); i++;
			label.replace(10*i, 2, "VY"); i++; break;
		case TVector3: label.replace(10*i, 2, "VX"); i++;
			label.replace(10*i, 2, "VY"); i++;
			label.replace(10*i, 2, "VZ"); i++; break;
		case TView: label.replace(10*i, 2, "VX"); i++;
			label.replace(10*i, 2, "VY"); i++;
			label.replace(10*i, 2, "VZ"); i++;
			label.replace(10*i, 3, "Ang"); i++; break;
		case TRGB: label.replace(10*i, 3, "Red"); i++;
			label.replace(10*i, 3, "Grn"); i++;
			label.replace(10*i, 3, "Blu"); i++; break;
		case TRGBA: label.replace(10*i, 3, "Red"); i++;
			label.replace(10*i, 3, "Grn"); i++;
			label.replace(10*i, 3, "Blu"); i++;
			label.replace(10*i, 3, "Alf"); i++; break;
		case TCMYK: label.replace(10*i, 3, "Cyn"); i++;
			label.replace(10*i, 3, "Mag"); i++;
			label.replace(10*i, 3, "Yel"); i++;
			label.replace(10*i, 3, "blK"); i++; break;
		default: label.replace(10*i, 4, "Real"); i++; break;
	}

	if ( p->next ) {
		order[7] = i;
		label.replace(10*i, 3, "FOM");
		i++;
	}
	label.replace(10*i, 1, "\0");

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG writeASCII: label: " << label << endl;
		cout << "DEBUG writeASCII: order: " << order[0] << " " << order[1] << " " << order[2] << " " << 
			order[3] << " " << order[4] << " " << order[5] << " " << order[6] << " " << order[7] << " " << order[8] << " (m=" << i << ")" << endl;
	}
	
	fimg << label << endl;
	
	long			c, x, y, z, n;
	float*			fom = (float *) p->next->data_pointer();
	
	for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			for ( y=0; y<p->sizeY(); y++ ) {
				for ( x=0; x<p->sizeX(); x++ ) {
					if ( order[0] > -1 ) fimg << " " << x;
					if ( order[1] > -1 ) fimg << " " << y;
					if ( order[2] > -1 ) fimg << " " << z;
					if ( order[3] > -1 ) fimg << " " << n;
					for ( c=0; c<p->channels(); c++, i++ ) fimg << " " << (*p)[i];
					if ( order[4] > -1 ) fimg << " " << fom[i];
					fimg << endl;
				}
			}
		}
	}
	
	fimg.close();
	
	return 0;
}

