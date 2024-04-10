/*
*
* Solver: QLPReader.cpp -- Copyright (c) 2010-2017 Jan Wolf
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "Utilities/ExternCode/QLPReader.hpp"

using namespace std;


string QLPReader::strtolower( const string& s ){
	string sl = s;
	transform( sl.begin(), sl.end(), sl.begin(), ::tolower );
	return sl;
}

// von http://www.c-plusplus.de/forum/125204-10
string QLPReader::trim( const string& s ){
	string::size_type first = s.find_first_not_of( " \n\t\r" );
	if( first == string::npos ) {
		return string();
	} else {
		string::size_type last = s.find_last_not_of( " \n\t\r" ); // must succeed
		return s.substr( first, last - first + 1 );
	}
}

string QLPReader::removeComments( const string& s ){
	string::size_type first = s.find_first_of( "\\" );
	if( first == string::npos ) {
		return s;
	} else {
		return s.substr( 0, first );
	}
}

const bool QLPReader::isNumeric( const std::string& s ){
	// TODO Check verbessern?

	const string tmp = trim( s );
	if( tmp.size() == 0 ){
		throw runtime_error( "Invalid isNumeric check." );
	}

	if( strtolower( s.substr( 0, 3 ) ) == "inf" ){
		return true;
	}

	const char t = s.at( 0 );
	return t== '+' || t == '-' || t == '0' || t == '1' || t == '2' || t == '3' || t == '4' || t == '5' || t == '6' || t == '7' || t == '8' || t == '9';
}



void QLPReader::readNumberAndShrink( stringHelper& s, mpq_class& n ){
	s.trim();

	// Vorzeichen
	bool signPos = true;

	while( s.length() ){
		// Schleife, um +- und -+ abzufangen
		if( s[0] == '+' ){
			s.removeFromBeg( 1 );
			s.trimLeft();
		} else if( s[0] == '-' ){
			signPos = false;
			s.removeFromBeg( 1 );
			s.trimLeft();
		} else {
			break;
		}
	}

	// Nur +/-, keine Zahl
	if( s.length() && !isdigit( s[0] ) ){
		if( !signPos ){
			n = mpq_class( -1, 1 );
		} else {
			n = mpq_class( 1, 1 );
		}
		return;
	}

	// (ersten Teil der) Zahl einlesen
	mpz_class n1 = 0;
	while( s.length() && isdigit( s[0] ) ){
		n1 *= 10;
		char t[2];
		t[0] = s[0];
		t[1] = 0;
		n1 += atoi( t );
		s.removeFromBeg( 1 );
	}

	// Je nach Zahl-Typ weitermachen
	if( s.length() > 1 && s[0] == '/' && ( isdigit( s[1] ) || s[1] == '+' || s[1] == '-' ) ){
		// Bruch
		s.removeFromBeg( 1 );

		// Vorzeichen Nenner
		bool sign2Pos = true;
		if( s[0] == '+' ){
			s.removeFromBeg( 1 );
			s.trimLeft();
		} else if( s[0] == '-' ){
			sign2Pos = false;
			s.removeFromBeg( 1 );
			s.trimLeft();
		}

		mpz_class n2 = 0;
		while( s.length() && isdigit( s[0] ) ){
			n2 *= 10;
			char t[2];
			t[0] = s[0];
			t[1] = 0;
			n2 += atoi( t );
			s.removeFromBeg( 1 );
		}
		if( signPos != sign2Pos ){
			n1 = -n1;
		}
		n = mpq_class( n1, n2 );
	} else if( s.length() && s[0] == '.' ){
		s.removeFromBeg( 1 );
		mpz_class n2 = 1;
		while( s.length() && isdigit( s[0] ) ){
			n1 *= 10;
			n2 *= 10;
			char t[2];
			t[0] = s[0];
			t[1] = 0;
			n1 += atoi( t );
			s.removeFromBeg( 1 );
		}
		if( !signPos ){
			n1 = -n1;
		}
		n = mpq_class( n1, n2 );
	} else {
		if( !signPos ){
			n1 = -n1;
		}
		n = n1;
	}

	// wissenschaftliche Darstellung
	if( s.length() > 1 && ( s[0] == 'e' || s[0] == 'E' ) && ( isdigit( s[1] ) || s[1] == '+' || s[1] == '-' ) ){
		s.removeFromBeg( 1 );
		bool sign3Pos = true;
		if( s[0] == '+' ){
			s.removeFromBeg( 1 );
		} else if( s[0] == '-' ){
			sign3Pos = false;
			s.removeFromBeg( 1 );
		}

		mpz_class n3 = 0;
		while( s.length() && isdigit( s[0] ) ){
			n3 *= 10;
			char t[2];
			t[0] = s[0];
			t[1] = 0;
			n3 += atoi( t );
			s.removeFromBeg( 1 );
		}

		// Schnelle Exponentation
		mpz_class n4 = 1;
		mpz_class x = 10;
		while( n3 > 0 ){
			if( n3 % 2 == 0 ){
				x *= x;
				n3 /= 2;
			} else {
				n4 *= x;
				n3 -= 1;
			}
		}

		if( !sign3Pos ){
			n /= n4;
		} else {
			n *= n4;
		}

	}

	n.canonicalize();
}


void QLPReader::readVariableAndShrink( stringHelper& s, const map<string, int>& varNumbers, int& n ){
	s.trim();
	string strn = s.substringUntilFirstOfAndShrink( " \n\t\r+-" );
	strn = trim( strn );

	auto it = varNumbers.find( strn );

	if( it == varNumbers.end() ){
		// Unexpected Variable
		throw string( strn );
	} else {
		n = it->second;
	}
}


void QLPReader::parseConstraint( const string& s, const map<string, int>& varNumbers, vector<constraintHelper>& matrix ){

	stringHelper tmpsh( s );
	tmpsh.trim();

	constraintHelper ch;

	tmpsh.removeUntilFirstOrDoNothing( ':' );
	// Namen einlesen und speichern?

	string s1 = tmpsh.substringUntilFirstOfAndShrink( "<>=" );

	if( tmpsh[0] == '=' ){
		tmpsh.removeFromBeg( 1 );
	}

	const char type = s1[ s1.length() - 1 ];
	stringHelper sh1( s1 );
	sh1.removeFromEnd( 1 );
	sh1.trim();
	tmpsh.trim();

	if( type == '<' ) ch.type = -1;
	else if( type == '>' ) ch.type = 1;
	else ch.type = 0;

	while( sh1.length() > 0 ){
		readerHelper h;
		readNumberAndShrink( sh1, h.mult );
		sh1.trimLeft();
		readVariableAndShrink( sh1, varNumbers, h.ind );
		sh1.trimLeft();
		if( h.mult != 0 ){
			ch.constraint.push_back( h );
		}
	}

	readNumberAndShrink( tmpsh, ch.rhs );
	matrix.push_back( ch );
}





data::Qlp QLPReader::read( const char* filename, const string randPrefix, const string allPrefix ){

	int objectiveBlockStart = -1;
	int objectiveBlockEnd = -1;

	int stBlockStart = -1;
	int stBlockEnd = -1;

	int boundBlockStart = -1;
	int boundBlockEnd = -1;

	int genBlockStart = -1;
	int genBlockEnd = -1;

	int binBlockStart = -1;
	int binBlockEnd = -1;

	int existsBlockStart = -1;
	int existsBlockEnd = -1;

	int allBlockStart = -1;
	int allBlockEnd = -1;

	int randomBlockStart = -1;
	int randomBlockEnd = -1;

	int ordBlockStart = -1;
	int ordBlockEnd = -1;

	map<string, int> varNumbers;
	vector<string> vars;
	vector<mpq_class> lbounds;
	vector<mpq_class> ubounds;
	vector<bool> linf;
	vector<bool> uinf;
	vector<string> specialBounds;
	vector<char> quantifiers;
	vector<char> numbersystems;

	bool maximize = false;
	vector<readerHelper> objfunc;
	vector<constraintHelper> matrix;


	vector<string> filecontent;


	// Datei einlesen
	{
		ifstream file( filename );
		string buffer = "";
		string tmpbuf = "";

		int* blockend = NULL;
		while( file.good() ){
			getline( file, buffer );
			buffer = removeComments( buffer );
			buffer = trim( buffer );
			if( buffer.empty() ){
				continue;
			}

			tmpbuf = strtolower( buffer.substr( 0, 10 ) );

			if( objectiveBlockStart == -1 && ( tmpbuf.substr( 0, 3 ) == "min" || tmpbuf.substr( 0, 3 ) == "max" ) ){
				objectiveBlockStart = filecontent.size();
				blockend = &objectiveBlockEnd;
			} else if( stBlockStart == -1 && ( tmpbuf.substr( 0, 10 ) == "subject to" || tmpbuf.substr( 0, 9 ) == "such that" || tmpbuf == "s.t." || tmpbuf == "st." || tmpbuf == "st" ) ){
				stBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &stBlockEnd;
			} else if( boundBlockEnd == -1 && ( tmpbuf == "bounds" || tmpbuf == "bound" ) ){
				boundBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &boundBlockEnd;
			} else if( genBlockEnd == -1 && ( tmpbuf == "generals" || tmpbuf == "general" || tmpbuf == "gen" ) ){
				genBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &genBlockEnd;
			} else if( binBlockEnd == -1 && ( tmpbuf == "binaries" || tmpbuf == "binary" || tmpbuf == "bin" ) ){
				binBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &binBlockEnd;
			} else if( existsBlockEnd == -1 && ( tmpbuf == "exists" || tmpbuf == "exist" ) ){
				existsBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &existsBlockEnd;
			} else if( allBlockEnd == -1 && ( tmpbuf == "all" ) ){
				allBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &allBlockEnd;
			} else if( randomBlockEnd == -1 && ( tmpbuf == "random" || tmpbuf == "rand" ) ){
				randomBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &randomBlockEnd;
			} else if( ordBlockEnd == -1 && ( tmpbuf == "order" || tmpbuf == "ord" ) ){
				ordBlockStart = filecontent.size();
				if( blockend ) *blockend = filecontent.size();
				blockend = &ordBlockEnd;
			} else if( tmpbuf == "end" ){
				if( blockend ) *blockend = filecontent.size();
			}

			filecontent.push_back( trim( buffer ) );

		}
	}



	// Bounds preprocessen
	if( boundBlockStart == -1 ){
		throw runtime_error( "Bounds missing." );
	}


	{
		unsigned int approxNumVars = boundBlockEnd - boundBlockStart + 1;

		vector<string> bndl;
		bndl.reserve( approxNumVars );
		vector<string> bndu;
		bndu.reserve( approxNumVars );
		vector<string> bndvar;
		bndvar.reserve( approxNumVars );
		vector<string> bndspec;
		bndspec.reserve( approxNumVars );

		for( int i = boundBlockStart + 1; i < boundBlockEnd; ++i ){
			string tmps = filecontent.at( i );

			string::size_type pos;
			while( ( pos = tmps.find( "<=" ) ) != string::npos ){
				tmps.erase( pos + 1, 1 );
			}
			while( ( pos = tmps.find( ">=" ) ) != string::npos ){
				tmps.erase( pos + 1, 1 );
			}
			while( ( pos = tmps.find( "==" ) ) != string::npos ){
				tmps.erase( pos + 1, 1 );
			}


			string::size_type eqpos = tmps.find_last_of( "=" );

			// Random-Spezialfall
			string::size_type spec = tmps.find_first_of( "[" );
			if( spec != string::npos ){
				if( eqpos == string::npos ){
					cerr << "Unable to parse: " << tmps << endl;
					throw runtime_error( "Unable to parse." );
				}

				string left = trim( tmps.substr( 0, eqpos ) );
				tmps.erase( 0, eqpos + 1 );
				string right = trim( tmps );

				bndl.push_back( "0" );
				bndvar.push_back( left );
				bndu.push_back( "0" );

				bndspec.resize( bndvar.size() );
				bndspec.at( bndvar.size() - 1 ) = right;
				continue;
			}


			string::size_type leqpos1 = tmps.find_first_of( "<" );
			string::size_type leqpos2 = tmps.find_last_of( "<" );

			string::size_type geqpos1 = tmps.find_first_of( ">" );
			string::size_type geqpos2 = tmps.find_last_of( ">" );

			if( leqpos1 != string::npos && leqpos1 != leqpos2 ){
				// l < x < u

				bndl.push_back( trim( tmps.substr( 0, leqpos1 ) ) );
				tmps.erase( 0, leqpos1 + 1 );
				leqpos2 = tmps.find_first_of( "<" );
				bndvar.push_back( trim( tmps.substr( 0, leqpos2 ) ) );
				tmps.erase( 0, leqpos2 + 1 );
				bndu.push_back( trim( tmps ) );

			} else if( geqpos1 != string::npos && geqpos1 != geqpos2 ){
				// u > x > l

				bndu.push_back( trim( tmps.substr( 0, geqpos1 ) ) );
				tmps.erase( 0, geqpos1 + 1 );
				geqpos2 = tmps.find_first_of( ">" );
				bndvar.push_back( trim( tmps.substr( 0, geqpos2 ) ) );
				tmps.erase( 0, geqpos2 + 1 );
				bndl.push_back( trim( tmps ) );

			} else if( leqpos1 != string::npos && leqpos1 == leqpos2 ){
				// l < x || x < u
				string l = trim( tmps.substr( 0, leqpos1 ) );
				tmps.erase( 0, leqpos1 + 1 );
				string u = trim( tmps );

				if( isNumeric( l ) && isNumeric( u ) ){
					cerr << "Unable to parse: " << u << " < " << right << endl;
					throw runtime_error( "Unable to parse." );
				}

				if( isNumeric( l ) ){
					bndl.push_back( l );
					bndvar.push_back( u );
					bndu.push_back( "inf" );
				} else if( isNumeric( u ) ){
					bndl.push_back( "0" );
					bndvar.push_back( l );
					bndu.push_back( u );
				} else {
					cerr << "Unable to parse: " << u << " < " << right << endl;
					throw runtime_error( "Unable to parse." );
				}

			} else if( geqpos1 != string::npos && geqpos1 == geqpos2 ){
				// l > x || x > u
				string u = trim( tmps.substr( 0, geqpos1 ) );
				tmps.erase( 0, geqpos1 + 1 );
				string l = trim( tmps );

				if( isNumeric( l ) && isNumeric( u ) ){
					cerr << "Unable to parse: " << u << " < " << right << endl;
					throw runtime_error( "Unable to parse." );
				}

				if( isNumeric( l ) ){
					bndl.push_back( l );
					bndvar.push_back( u );
					bndu.push_back( "inf" );
				} else if( isNumeric( u ) ){
					bndl.push_back( "0" );
					bndvar.push_back( l );
					bndu.push_back( u );
				} else {
					cerr << "Unable to parse: " << u << " < " << right << endl;
					throw runtime_error( "Unable to parse." );
				}

			} else if( eqpos != string::npos ){
				// x = const

				string left = trim( tmps.substr( 0, eqpos ) );
				tmps.erase( 0, eqpos + 1 );
				string right = trim( tmps );

				if( isNumeric( left ) && !isNumeric( right ) ){
					bndl.push_back( left );
					bndvar.push_back( right );
					bndu.push_back( left );
				} else if( isNumeric( right ) && !isNumeric( left ) ){
					bndl.push_back( right );
					bndvar.push_back( left );
					bndu.push_back( right );
				} else {
					cerr << "Unable to parse: " << left << " = " << right << endl;
					throw runtime_error( "Unable to parse." );
				}

			} else {

				if( tmps.find_first_of( " " ) == tmps.find_last_of( " " ) && tmps.find_first_of( " " ) != string::npos ){
					// Free variable
					bndl.push_back( "-inf" );
					bndvar.push_back( tmps.substr( 0, tmps.find_first_of( " " ) ) );
					bndu.push_back( "inf" );
				} else {
					cerr << tmps << endl;
					throw runtime_error( "Unable to read bound." );
				}

			}
		}

		bndspec.resize( bndvar.size() );




		// Variablen finden
		{
			if( ordBlockStart != -1 ){
				// Variablenreihenfolge aus Order-Block auslesen
				for( int i = ordBlockStart + 1; i < ordBlockEnd; ++i ){
					stringstream ss( filecontent.at( i ) );
					string temp;
	//				while( ss && getline( ss, temp, ' ' ) ){
					while( ss && ( ss >> temp ) ){
						if( varNumbers.find( temp ) != varNumbers.end() ){
							throw runtime_error( "Variable doppelt in ORDER-Block: " + temp );
						}
						varNumbers[ temp ] = vars.size();
						vars.push_back( temp );
					}
				}
			} else {
				// Variablenreihenfolge aus Bounds-Block auslesen
				for( unsigned int i = 0; i < bndvar.size(); ++i ){
					const string temp = bndvar.at( i );
					if( varNumbers.find( temp ) != varNumbers.end() ){
						throw runtime_error( "Variable doppelt in ORDER-Block." );
					}
					varNumbers[ temp ] = vars.size();
					vars.push_back( temp );
				}
			}
		}





		// Bounds lesen
		{
			lbounds.resize( vars.size() );
			ubounds.resize( vars.size() );
			linf.resize( vars.size(), false );
			uinf.resize( vars.size(), false );
			specialBounds.resize( vars.size() );

			if( bndvar.size() != vars.size() ){
				throw runtime_error( "Wrong number of variables." );
			}

			for( unsigned int i = 0; i < vars.size(); ++i ){

				if( varNumbers.find( bndvar.at( i ) ) == varNumbers.end() ){
					throw runtime_error( "Error in variables." );
				}

				mpq_class b;

				stringHelper lbndsh( bndl.at( i ) );
				readNumberAndShrink( lbndsh, b );
				lbndsh.trim();

				if( lbndsh.length() ){
					linf.at( varNumbers.at( bndvar.at( i ) ) ) = true;
	//				cout << bndl.at( i ) << endl;
	//				throw runtime_error( "Error in lower bounds." );
				} else {
					lbounds.at( varNumbers.at( bndvar.at( i ) ) ) = b;
				}

				stringHelper rbndsh( bndu.at( i ) );
				readNumberAndShrink( rbndsh, b );
				rbndsh.trim();

				if( rbndsh.length() ){
					uinf.at( varNumbers.at( bndvar.at( i ) ) ) = true;
	//				cout << bndu.at( i ) << endl;
	//				throw runtime_error( "Error in upper bounds." );
				} else {
					ubounds.at( varNumbers.at( bndvar.at( i ) ) ) = b;
				}

				specialBounds.at( varNumbers.at( bndvar.at( i ) ) ) = bndspec.at( i );
			}
		}

	}




	// Objective einlesen
	{
		if( strtolower( filecontent.at( objectiveBlockStart ).substr( 0, 3 ) ) == "max"){
			maximize = true;
		}

		string objectiveStr = "";
		for( int i = objectiveBlockStart + 1; i < objectiveBlockEnd; ++i ){
			objectiveStr += filecontent.at( i );
			objectiveStr += " ";
		}

		stringHelper objective( objectiveStr );
		objective.removeUntilFirstOrDoNothing( ':' );

		objective.trim();

		while( objective.length() > 0 ){
			readerHelper h;
			readNumberAndShrink( objective, h.mult );
			try{
				readVariableAndShrink( objective, varNumbers, h.ind );
			} catch( string &var ){
				lbounds.push_back( 0 );
				ubounds.push_back( 0 );
				linf.push_back( false );
				uinf.push_back( true );
				h.ind = vars.size();
				varNumbers[ var ] = vars.size();
				vars.push_back( var );
				cerr << "Hinweis: Vorher unbekannte Variable in der Zielfunktion: " << var << endl;
			}
			objfunc.push_back( h );
		}
	}





	// Constraints einlesen
	if( stBlockStart != -1 ){

		string constraint;

		for( int i = stBlockStart + 1; i < stBlockEnd; ++i ){
			constraint += filecontent.at( i );
			constraint += " ";

			string::size_type tmppos = filecontent.at( i ).find_last_of( "<>=" );
			if( tmppos != string::npos ) {
				while( !constraint.empty() ){
					try{
						parseConstraint( constraint, varNumbers, matrix );
						constraint.clear();
					} catch( string &var ){
						lbounds.push_back( 0 );
						ubounds.push_back( 0 );
						linf.push_back( false );
						uinf.push_back( true );
						varNumbers[ var ] = vars.size();
						vars.push_back( var );
						cerr << "Hinweis: Vorher unbekannte Variable in den Nebenbedingungen: " << var << endl;
					}
				}
			}
		}
	}




	// Quantoren lesen
	{
		quantifiers.resize( vars.size() );
		for( unsigned int i = 0; i < vars.size(); ++i ){
			quantifiers.at( i ) = 'E';
		}

		// Allquantoren
		if( allBlockStart != -1 ){
			string alls = "";
			for( int i = allBlockStart + 1; i < allBlockEnd; ++i ){
				alls += filecontent.at( i );
				alls += " ";
			}
			stringstream allss( alls );
			try {
				for( std::string var; allss >> var; quantifiers.at( varNumbers.at( var ) ) = 'A' );
			} catch( out_of_range &e ){
				throw runtime_error( "Invalid All Variable." );
			}

		}

		// Randomquantoren
		if( randomBlockStart != -1 ){
			string rands = "";
			for( int i = randomBlockStart + 1; i < randomBlockEnd; ++i ){
				rands += filecontent.at( i );
				rands += " ";
			}
			stringstream randss( rands );
			try {
				for( std::string var; randss >> var; quantifiers.at( varNumbers.at( trim( var ) ) ) = 'R' );
			} catch( out_of_range &e ){
				throw runtime_error( "Invalid Random Variable." );
			}
		}

		// Ggf. Prefixe verarbeiten

		if( existsBlockStart == -1 && allBlockStart == -1 && randomBlockStart == -1 && !randPrefix.empty() ){
			for( unsigned int i = 0; i < vars.size(); ++i ){
				if( vars.at( i ).substr( 0, randPrefix.length() ) == randPrefix ){
					quantifiers.at( i ) = 'R';
				}
			}
		}

		if( existsBlockStart == -1 && allBlockStart == -1 && randomBlockStart == -1 && !allPrefix.empty() ){
			for( unsigned int i = 0; i < vars.size(); ++i ){
				if( vars.at( i ).substr( 0, allPrefix.length() ) == allPrefix ){
					quantifiers.at( i ) = 'A';
				}
			}
		}
	}



	// Numbersystem lesen
	{
		numbersystems.resize( vars.size() );
		for( unsigned int i = 0; i < vars.size(); ++i ){
			numbersystems.at( i ) = 'R';
		}

		// Generals
		if( genBlockStart != -1 ){
			string gens = "";
			for( int i = genBlockStart + 1; i < genBlockEnd; ++i ){
				gens += filecontent.at( i );
				gens += " ";
			}
			stringstream genss( gens );
			try {
				for( std::string var; genss >> var; numbersystems.at( varNumbers.at( var ) ) = 'G' );
			} catch( out_of_range &e ){
				throw runtime_error( "Invalid Integer Variable." );
			}
		}

		// Binaries
		if( binBlockStart != -1 ){
			string bins = "";
			for( int i = binBlockStart + 1; i < binBlockEnd; ++i ){
				bins += filecontent.at( i );
				bins += " ";
			}
			stringstream binss( bins );
			try {
				for( std::string var; binss >> var; numbersystems.at( varNumbers.at( var ) ) = 'B' );
			} catch( out_of_range &e ){
				throw runtime_error( "Invalid Binary Variable." );
			}
		}
	}




	// in QLP-Objekt umwandeln

	vector<data::QpVar> qpvars;
	for( unsigned int i = 0; i < vars.size(); ++i ){
		data::QpVar::NumberSystem ns;

		switch( numbersystems.at( i ) ){
			case 'R':
				ns = data::QpVar::real;
				break;
			case 'G':
				ns = data::QpVar::generals;
				break;
			case 'B':
				ns = data::QpVar::binaries;
				break;
			default:
				throw runtime_error( "unexpected error" );
				break;
		}

		data::QpVar::Quantifier q;

		switch( quantifiers.at( i ) ){
			case 'E':
				q = data::QpVar::exists;
				break;
			case 'A':
				q = data::QpVar::all;
				break;
			case 'R':
				q = data::QpVar::random;
				break;
			default:
				throw runtime_error( "unexpected error" );
				break;
		}

		if( q == data::QpVar::random ){

			if( specialBounds.at( i ).empty() ){

				// Distributions nicht angegeben?
				// Gleichverteilung �ber bounds

				if( 1 ){
					throw runtime_error( "Bound-Error. No Dirstribution given. Error. Finish." ); exit(0);
				}

			} else {
				// Specialbounds einlesen

				string::size_type tmpstart = specialBounds.at( i ).find_first_of( "[" );
				string::size_type tmpend = specialBounds.at( i ).find_first_of( "]" );

				if( tmpstart == string::npos || tmpend == string::npos ) {
					throw runtime_error( "Bound-Range-Error." );
				}

				string rangestring = trim( specialBounds.at( i ).substr( tmpstart + 1, tmpend - tmpstart - 1 ) );


				string distrstring = specialBounds.at( i ).substr( tmpend + 1, string::npos );

				tmpstart = distrstring.find_first_of( "[" );
				tmpend = distrstring.find_first_of( "]" );

				if( tmpstart == string::npos || tmpend == string::npos ) {
					throw runtime_error( "Bound-Range-Error." );
				}

				distrstring = trim( distrstring.substr( tmpstart + 1, tmpend - tmpstart - 1 ) );


				mpq_class tmpn;

				vector<data::QpNum> rangevec;
				stringHelper rngStrH( rangestring );
				while( rngStrH.length() ){
					readNumberAndShrink( rngStrH, tmpn );
					rangevec.push_back( (double)tmpn.getNominator() / (double)tmpn.getDenominator() );
					rngStrH.trim();
					while( rngStrH.length() && rngStrH[0] == ',' ){
						rngStrH.removeFromBeg( 1 );
						rngStrH.trimLeft();
					}
				}

				vector<data::QpRational> distrvec;
				stringHelper distStrH( distrstring );
				while( !distrstring.empty() ){
					readNumberAndShrink( distStrH, tmpn );
					distrvec.push_back( data::QpRational(tmpn.getNominator(),tmpn.getDenominator()) );
					distStrH.trim();
					while( distStrH.length() && distStrH[0] == ',' ){
						distStrH.removeFromBeg( 1 );
						distStrH.trimLeft();
					}
				}

				qpvars.push_back( data::QpVar( vars.at( i ), i, rangevec, distrvec, ns ) );
			}
		} else {
			data::QpNum l;
			data::QpNum u;

			if( linf.at( i ) ){
				l = data::QpNum( true );
			} else {
				l = (double)lbounds.at(i).getNominator() / (double)lbounds.at(i).getDenominator();
			}

			if( uinf.at( i ) ){
				u = data::QpNum( false );
			} else {
				u = (double)ubounds.at( i ).getNominator() / (double)ubounds.at( i ).getDenominator();
			}

			qpvars.push_back( data::QpVar( vars.at( i ), i, l, u, ns, q ) );
		}
	}


	data::QpObjFunc obj;
	if( maximize ){
		obj.setObjective( data::QpObjFunc::max );
	}
	obj.setSize( vars.size() );
	for( unsigned int i = 0; i < objfunc.size(); ++i ){
		obj.setObjElement( objfunc.at( i ).ind, data::QpNum( (double)objfunc.at( i ).mult.getNominator() / (double)objfunc.at( i ).mult.getDenominator() ) );
	}


	vector<data::QpRhs> rhs;

	for( unsigned int i = 0; i < matrix.size(); ++i ){
		data::QpRhs::RatioSign s;
		switch( matrix.at( i ).type ){
			case -1:
				s = data::QpRhs::smallerThanOrEqual;
				break;
			case 0:
				s = data::QpRhs::equal;
				break;
			case 1:
				s = data::QpRhs::greaterThanOrEqual;
				break;
			default:
				throw runtime_error( "unexpected error" );
				break;
		}
		rhs.push_back( data::QpRhs( data::QpNum( (double)matrix.at( i ).rhs.getNominator() / (double)matrix.at( i ).rhs.getDenominator()), s ) );
	}


//	vector<data::QpSparseVec> lhs( matrix.size(), vars.size() );
	data::QpSparseMatrix lhs( matrix.size() );

	for( unsigned int i = 0; i < matrix.size(); ++i ){
		for( unsigned int j = 0; j < matrix.at( i ).constraint.size(); ++j ){
			lhs.at( i ).push_back( data::IndexedElement( matrix.at( i ).constraint.at( j ).ind, data::QpNum( (double)matrix.at( i ).constraint.at( j ).mult.getNominator() / (double)matrix.at( i ).constraint.at( j ).mult.getDenominator() ) ) );
//			lhs.at( i )[matrix.at( i ).constraint.at( j ).ind] = data::QpNum( matrix.at( i ).constraint.at( j ).mult );
		}
	}

//	Qlp(const data::QpObjFunc& objective,
//			const std::vector<data::QpVar>& varVec,
//			const data::QpSparseMatrix&,
//			const std::vector<data::QpRhs>& rhsVec);

	data::Qlp qlp( obj, qpvars, lhs, rhs );

	return qlp;
}
