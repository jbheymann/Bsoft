		if ( curropt->tag == "Volt" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-Volt: A voltage must be specified!" << endl;
			else
				jsctf["Volt"] = v;
		}
		if ( curropt->tag == "Amplitude" ) {
			v = curropt->value.real();
			if ( v < -1 || v > 1 )
				cerr << "-Amplitude: A fraction must be specified!" << endl;
			else
				jsctf["Amplitude_phase"] = asin(v);
		}
		if ( curropt->tag == "Cs" ) {
			v = curropt->real_units();
			jsctf["Cs"] = v;
		}
		if ( curropt->tag == "Cc" ) {
			if ( ( v = curropt->real_units() ) < 0.1 )
				cerr << "-Cc: A Cc value must be specified!" << endl;
			else
				jsctf["Cc"] = v;
		}
		if ( curropt->tag == "alpha" ) {
			if ( ( v = curropt->real_units() ) < 1e-10 )
				cerr << "-alpha: A value must be specified!" << endl;
			else {
				jsctf["Beam_convergence"] = v;
			}
		}
		if ( curropt->tag == "energyspread" ) {
			if ( ( v = curropt->value.real() ) < 1e-10 )
				cerr << "-energyspread: A value must be specified!" << endl;
			else
				jsctf["Energy_spread"] = v;
		}
		if ( curropt->tag == "aperture" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-aperture: An aperture diameter must be specified!" << endl;
			else
				jsctf["Objective_aperture"] = v;
		}
		if ( curropt->tag == "focallength" ) {
			if ( ( v = curropt->real_units() ) < 1 )
				cerr << "-focallength: A focal length must be specified!" << endl;
			else
				jsctf["Focal_length"] = v;
		}
		if ( curropt->tag == "slitwidth" ) {
			if ( ( v = curropt->value.real() ) < 1 )
				cerr << "-slitwidth: A slit width must be specified!" << endl;
			else
				jsctf["Slit_width"] = v;
		}
		if ( curropt->tag == "basetype" ) {
			long	basetype = curropt->value.integer();
			if ( basetype < 1 || basetype > 6 ) {
				basetype = 1;
				cerr << "Warning: The baseline type must be 1-6. Reset to 1." << endl;
			} else
				jsctf["Baseline_type"] = basetype;
		}
		if ( curropt->tag == "baseline" ) {
			vector<double>	d = parse_real_vector(curropt->value.c_str());
			if ( d.size() < 1 )
				cerr << "-baseline: At least one coefficient must be specified!" << endl;
			else
				jsctf["Baseline"] = d;
		}
		if ( curropt->tag == "envtype" ) {
			long	envtype = curropt->value.integer();
			if ( envtype < 1 || envtype > 4 ) {
				envtype = 4;
				cerr << "Warning: The envelope type must be 1-4. Reset to 4." << endl;
			} else
				jsctf["Envelope_type"] = envtype;
		}
		if ( curropt->tag == "envelope" ) {
			vector<double>	d = parse_real_vector(curropt->value.c_str());
			if ( d.size() < 1 )
				cerr << "-envelope: At least an envelope amplitude must be specified!" << endl;
			else
				jsctf["Envelope"] = d;
		}
