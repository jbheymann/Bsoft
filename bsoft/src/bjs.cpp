/**
@file	bjs.cpp
@author	Bernard Heymann
@date	Created: 20160911
@date	Modified: 20240728

clang++ -o js *.cpp -std=c++11

**/

#include "json.h"
#include "STARparser.h"
#include "mdoc.h"
#include "utilities.h"


template <typename T>
ostream&	operator<<(ostream& output, vector<T> vec) {
	for ( auto it = vec.begin(); it != vec.end(); ++it )
		cout << "\t" << *it;
	return output;
}

void		show_type_sizes()
{
	cout << "JStype size: " << sizeof(JStype) << endl;
	cout << "string size: " << sizeof(string) << endl;
	cout << "bool size: " << sizeof(bool) << endl;
	cout << "long size: " << sizeof(long) << endl;
	cout << "double size: " << sizeof(double) << endl;
	cout << "vector<JSvalue> size: " << sizeof(vector<JSvalue>) << endl;
	cout << "map<string, JSvalue> size: " << sizeof(map<string, JSvalue>) << endl;
	cout << "JSvalue size: " << sizeof(JSvalue) << endl;
}

JSvalue		readJSON(string filename)
{
	string			ext(extension(filename));

	JSvalue			root;
	
	if ( ext == "json" ) {
		JSparser		parser(filename);
		root = parser.parse();
	} else if ( ext == "star" ) {
		STARparser		parser(filename);
		root = parser.parse();
	} else if ( ext == "mdoc" ) {
		MDOCparser		parser(filename);
		root = parser.parse();
	} else {
		cerr << "Error: Extension " << ext << " not supported!" << endl;
	}
	
	return root;
}

int			writeJSON(string filename, JSvalue& root)
{
	string	ext = extension(filename);
	cout << "Writing " << filename << endl;
	
	if ( ext == "json" ) {
		ofstream        f(filename);
    	if ( f.fail() ) {
			cerr << "File " << filename << " not opened!" << endl;
			return -1;
		}
		
		f << root << endl;
//		f << ss.rdbuf();
		
		f.close();
		
	} else if ( ext == "star" ) {
		writeSTAR(filename, root);
	} else {
		cerr << "Error: Extension " << ext << " not supported!" << endl;
	}

	return 0;
}

int		mg_convert(JSvalue& root)
{
	JSvalue				parr(JSarray);
	
	for ( auto it = root["data"].begin(); it != root["data"].end(); ++it ) {
		JSvalue&		mgfn = (*it)["micrograph.file_name"];
		if ( mgfn.type() == JSnull ) mgfn = (*it)["particle.file_name"];
		cout << "mg: " << mgfn << " (" << mgfn.type() << ")" << endl;
		for ( auto ip: (*it)["particle"] ) {
			ip["micrograph.file_name"] = mgfn;
			parr.push_back(ip);
		}
		it->erase("particle");
	}

	root["particle"] = parr;

	root["micrograph"] = root["data"];
	
	root.erase("data");
	
	return 0;
}

int		mod_convert(JSvalue& root)
{
	root["model"] = root["data"];
	
	root.erase("data");
	
	return 0;
}

// Usage assistance
const char* use[] = {
" ",
"Usage: bjs [options] in.star out.json",
"-------------------------------------",
"For JSON files.",
" ",
"Actions:",
"-show                    Print contents.",
"-convert mg              Convert a micrograph or model.",
" ",
NULL
};

/*
	JSON path query syntax:
		.
		"[]"
		"$['JSONimage']['Size'][0]"
		"['JSONimage']['Size'][0]"
		"$.JSONimage.Size[1]"
		".JSONimage.Size[1]"
		.JSONimage.Size

*/
int		main(int argc, char** argv)
{
	if ( argc < 2 ) {
		for ( long i = 0; use[i] != NULL; ++i )
			cout << use[i] << endl;
		return 0;
	}
	
	bool			show(0);
	long			optind, setstr(0), setint(0), setreal(0), integer(0);
	string			convert;
	double			real(0);
	vector<long>	intvec;
	vector<double>	realvec;
	string			jsonpath, str;
	
	for ( optind=1; optind<argc && argv[optind][0] == '-'; ++optind ) {
//		cout << argv[optind] << endl;
		if ( strcmp(argv[optind], "-show") == 0 ) show = 1;
		if ( strcmp(argv[optind], "-types") == 0 ) show_type_sizes();
		if ( strcmp(argv[optind], "-path") == 0 ) jsonpath = argv[++optind];
		if ( strcmp(argv[optind], "-set") == 0 ) {
			setstr = 1;
			str = argv[++optind];
		}
		if ( strcmp(argv[optind], "-setint") == 0 ) {
			setint = 1;
			integer = stol(argv[++optind]);
		}
		if ( strcmp(argv[optind], "-setreal") == 0 ) {
			setreal = 1;
			real = stod(argv[++optind]);
		}
		if ( strcmp(argv[optind], "-setintvec") == 0 ) intvec = parse_integer_vector(argv[++optind]);
		if ( strcmp(argv[optind], "-setrealvec") == 0 ) realvec = parse_real_vector(argv[++optind]);
		if ( strcmp(argv[optind], "-convert") == 0 )
			convert = argv[++optind];
	}
	
	string			filename(argv[optind++]);

	JSvalue			root = readJSON(filename);
		
	JSvalue			rootcopy = root;
	
	if ( jsonpath.length() ) {
		vector<JSvalue*>	valarr = root(jsonpath);
		if ( setstr ) *valarr[0] = str;
		if ( setint ) *valarr[0] = integer;
		if ( setreal ) *valarr[0] = real;
		if ( intvec.size() ) *valarr[0] = JSvalue(intvec);
		if ( realvec.size() ) *valarr[0] = JSvalue(realvec).array();
		cout << jsonpath << ": " << *valarr[0] << endl;
		cout << jsonpath << ": " << *rootcopy(jsonpath)[0] << endl;
	}
	
	if ( convert == "mg" )
		mg_convert(root);
	else if ( convert == "mod" )
		mod_convert(root);
	
	if ( show ) cout << root << endl;
	
	cout << "Memory: " << root.memory() << endl;
	
	stringstream	ss;
	ss << root << endl;
	cout << "Size: " << ss.tellp() << endl;

	if ( optind < argc ) {
		filename = argv[optind];
		writeJSON(filename, root);
	}
	
	return 0;
}

