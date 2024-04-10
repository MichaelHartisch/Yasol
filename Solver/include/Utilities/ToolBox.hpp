/*
*
* Solver: Toolbox.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef TOOLBOX_HPP_
#define TOOLBOX_HPP_

#include "Settings/Settings.hpp"

#include <map>
#include <iostream>
#include <cmath>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

//#include <gmpxx.h>
//#include <boost/filesystem.hpp>
//#include <boost/filesystem/fstream.hpp>
#include <map>
#include <algorithm>
#include <list>
#include <vector>

//#include "TOSimplex/TOFileReaderLP.h"
#include "Datastructures/numbers/myRatNum.hpp"


#include <cstdlib>
#include<iostream>
#include<fstream>

namespace utils {
class ToolBox {

public:

//	typedef enum {
//		IN, OUT
//	} FileType;
//
    
	static void readFromFile(const std::string& path, std::list<std::string>& target) {
		try {
			/*
			boost::filesystem::path p(path);
			if (boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p)) {
				boost::filesystem::ifstream ifs(p);
				target.clear();
				std::string line;
				while (!ifs.eof()) {
					getline(ifs, line);
					if (line.size() > 0) {
						target.push_back(line);
					}
				}
				ifs.close();
				*/
	        std::ifstream ifs(path.c_str());
	        std::string line;
	        //for (string s;!filestream.eof();) w�re wohl eine Alternative, soll aber nicht so gut sein.
	        while(ifs.good()) {
	            getline(ifs,line);
				if (line.size() > 0) {
					target.push_back(line);
				}
	        }
	        ifs.close();
		}
		catch (std::ios_base::failure &fail ) {
			//std::cerr << e.what() << std::endl;
		}
	}

	static std::list<std::string> readFromFile(const std::string& path) {
		std::list<std::string> stringList;
		utils::ToolBox::readFromFile(path, stringList);
		return stringList;
	}

	static void writeToFile(const std::string& filename, const std::string& content) {
		try {
			std::ofstream ofs;
			ofs.open(filename.c_str());
			ofs << content;
			ofs.close();
		}
		catch (std::ios_base::failure &fail) {
			//std::cerr << e.what() << std::endl;
		}
	}

	static void getDirContent(std::string dir, std::vector<std::string> &files) {
		files.clear();
		try {
		    DIR           *dirHandle;
		    struct dirent *dirEntry;

		    dirHandle = opendir(dir.c_str());
		    if ( dirHandle != NULL )
		    {
		        while ( 0 != ( dirEntry = readdir( dirHandle ) ) )
		        {
		            struct stat entrystat;
		            if (stat(dirEntry->d_name, &entrystat)) {
		                continue;
		            }
		        	if (S_ISREG(entrystat.st_mode)) {
					    files.push_back(std::string(dirEntry->d_name));
		        	}
		        }
		        closedir( dirHandle );
		    }
		}
		catch (std::ios_base::failure &fail) {
			//std::cerr << e.what() << std::endl;
		}
	}

	static std::map<std::string, std::string> readConfigFile(std::string path) {
		std::map<std::string, std::string> hashTable;
		std::list<std::string> stringList = readFromFile(path);
		std::list<std::string>::iterator start = stringList.begin();
		std::list<std::string>::iterator end = stringList.end();
		//For all lines in the file
		std::string::size_type pos;
		std::string s1;
		std::string s2;
		while (start != end) {
			//Parse Lines and throw away malformed lines: Entries consist of : key = value
			//For every line add key/value pair in map
			pos = (*start).find_first_of("=");
			if (pos != std::string::npos) {
				s1 = (*start).substr(0, pos);
				s2 = (*start).substr(pos + 1);
				hashTable.insert(make_pair(utils::ToolBox::removeStartEndWhitespaces(s1), utils::ToolBox::removeStartEndWhitespaces(s2)));
			}
			start++;
		}
		return hashTable;
	}

	//---------------------
	// von http://www.c-plusplus.de/forum/125204-10
	static const std::string trim(const std::string& s) {
		std::string::size_type first = s.find_first_not_of(" \n\t\r");
		if (first == std::string::npos) {
			return std::string();
		} else {
			std::string::size_type last = s.find_last_not_of(" \n\t\r"); // must succeed
			return s.substr(first, last - first + 1);
		}
	}

	/*static RatInf<mpq_class> stringToTORationalInf(const std::string& s){
		TOFileReaderLP r;
		RatInf<mpq_class> n;
		//std::cerr << "read number: " << s << std::endl;
		r.readInfNumber(s,n);
		return n;
	}*/

	static mpq_class stringToMpqClass(const std::string& s) {
		mpq_class value(0,0);
		mpq_class minus1(-1,1);
		std::string tmps = trim(s);

		bool signPos = true;
		if (tmps[0] == '+') {
			tmps.erase(0, 1);
			trim(s);
		} else if (tmps[0] == '-') {
			signPos = false;
			tmps.erase(0, 1);
			trim(s);
		}

		mpq_class n1 = 0, ten = 10;
		while (isdigit(tmps[0])) {
			n1 = n1 * ten;
			char t[2];
			t[0] = tmps[0];
			t[1] = 0;
			mpq_class add = atoi(t);
			n1 = n1 + add;
			tmps.erase(0, 1);
		}

		//std::cerr << "A extracted: " << n1.getNominator() << " / " << n1.getDenominator() << std::endl;
		//std::cerr << "B extract: " << tmps << std::endl;
		if (tmps[0] == '.') {
			tmps.erase(0, 1);
			mpq_class n2 = 1;
			while (isdigit(tmps[0])) {
				n1 = n1 * ten;
				n2 = n2 * ten;
				char t[2];
				t[0] = tmps[0];
				t[1] = 0;
				mpq_class add = atoi(t);
				n1 = n1 + add;
				tmps.erase(0, 1);
			}

			if (!signPos) {
				n1 = n1 * minus1;
			}
			value = n1 / n2;
		} else {
			if (!signPos) {
				n1 = n1 * minus1;
			}
			value = n1;
		}
		if (tmps.size() >0){
		    if (tmps[0] == 'e' || tmps[0] == 'E'){
                    tmps.erase(0, 1);

                    //Assume that such numbers are written without any spaces between e, +/- and the exponent
                signPos=true;
                if (tmps[0] == '-') signPos=false;
                if (tmps[0] == '-' || tmps[0] == '+') tmps.erase(0, 1);
                mpq_class exponent=0;
                while(tmps.size()>0){
                    if(!isdigit(tmps[0])){
                            std::cerr << "QpNum Warning: Exponential of Scientific QpNum is not integer: " << tmps << std::endl;
    //	                        assert(0);
                    }
                    exponent = exponent * ten;
	                char t[2];
                    t[0] = tmps[0];
                    t[1] = 0;
                    mpq_class add = atoi(t);
                    exponent = exponent + add;
                    tmps.erase(0, 1);
                }
                while(exponent>0){
                    if(signPos) value =value*ten;
                    else value = value/ten;
                    exponent = exponent+minus1;
                }
            }
		    else{
			std::cerr << "QpNum Warning: Remainder of QpNum is not exponential: " << tmps << std::endl;
//			assert(0)
		    }
		}
		//std::cerr << "B extracted: " << value.getNominator() << " / " << value.getDenominator() << std::endl;

		return value;
	}

	// -------------------------- Coversion of some Datastructures to a string representation --------------->
	static std::string boolToString(int b) {
		if (!b)
			return "false";
		return "true";
	}

//	static std::string mpq_classVecToString(const std::vector<mpq_class>& vec) {
//		std::string s("< ");
//		for (unsigned int i = 0; i < vec.size(); i++) {
//			s += convertToString(vec[i].get_d());
//			if (i != vec.size() - 1)
//				s += ", ";
//		}
//		s += " >";
//		return s;
//	}
//
//	static std::string mpq_classToString(const mpq_class& elem) {
//		return convertToString(elem.get_d());
//	}

	template<typename T> static std::string vecToString(const std::vector<T>& vec, int precision = TO_STRING_ACCURACY) {
		std::string begin("< ");
		unsigned int size = vec.size();
		for (unsigned int i = 0; i < size; i++) {
			begin += utils::ToolBox::convertToString(vec[i], precision);
			if (i < size - 1)
				begin += ", ";
		}
		begin += " >";
		return begin;
	}

	template<typename T> static std::string pairToString(const std::pair<T, T>& p) {
		std::string s("<");
		s += convertToString(p.first);
		s += "," + convertToString(p.second);
		s += ">";
		return s;
	}

	template<typename T>
	T static StringToNumber(const std::string &Text, T defValue = T()) {
		std::stringstream ss;
		for (std::string::const_iterator i = Text.begin(); i != Text.end(); ++i)
			if (isdigit(*i) || *i == 'e' || *i == '-' || *i == '+' || *i == '.')
				ss << *i;
		T result;
		if (ss >> result) {
			return result;
		} else {
			std::cout << "-----------------------------------------> StringToNumber conversion failed: " << Text << std::endl;
			return defValue;
		}
	}

	template<typename T> static std::string pairListToString(const std::list<std::pair<T, T> >& pList) {
		std::string s("[");
		std::list<std::pair<int, int> >::iterator it = pList.begin();
		while (it != pList.end()) {
			s += pairToString((*it));
			it++;
		}
		s += "]";
		return s;
	}

	template<typename T> static std::string pairVecToString(const std::vector<std::pair<T, T> >& vec) {
		std::string begin("< ");
		for (unsigned int i = 0; i < vec.size(); i++) {
			begin += utils::ToolBox::pairToString(vec[i]);
			if (i < vec.size() - 1)
				begin += ",";
		}
		begin += " >";
		return begin;
	}

	template<class T> static std::string to_string(const T& t) {
		std::stringstream ss;
		ss << t;
		return ss.str();
	}

	template<typename T> static std::string convertToString(const T& val, int precision = TO_STRING_ACCURACY) {
		if (!(val > MAX_QPDOUBLE_INF || val < MIN_QPDOUBLE_INF)) {
			if (val == MIN_QPDOUBLE_INF)
				return std::string("-inf");
			if (val == MAX_QPDOUBLE_INF)
				return std::string("+inf");
		}
		std::ostringstream out;
		T v = val;
		if((T)(int)v==val)
			precision=0;
		out.precision(precision);
		if(EXP_NOTATION){
			out << val;
		}else{
			out << std::fixed << val;
		}



		return (out.str());
	}

	static std::string listToString(const std::list<double>& list) {
		std::string s("< ");
		unsigned int i = 0;
		std::list<double>::const_iterator it = list.begin();
		while (it != list.end()) {
			s += convertToString(*it);
			if (i < list.size() - 1)
				s += ",";
			it++;
			i++;
		}
		s += " >";
		return s;
	}

	static std::string& trimString(std::string& s) {
		removeComments(s);
		return removeStartEndWhitespaces(s);
	}

	static std::string& removeStartEndWhitespaces(std::string& s) {
		unsigned int i = s.size();
		if (!i)
			return s;
		unsigned int j = 0;
		// check for white spaces at the beginning of the line
		while (s[j] == ' ')
			j++;
		// check for white spaces at the end of the line
		while (i > 0 && s[i - 1] == ' ')
			i--;
		// remove white spaces at the beginning and the end of the line
		if (i < s.size())
			s = s.substr(0, i);
		if (j > 0 && s.size())
			s = s.substr(j);
		return s;
	}

	static std::string& removeWhitespaces(std::string& str) {
		for (unsigned int i = 0; i < str.length(); i++)
			if (str[i] == ' ') {
				str.erase(i, 1);
				i--;
			}
		return str;
	}

	static std::string& removeComments(std::string& s) {
		std::string::size_type pos0;
		if ((pos0 = s.find("//")) != std::string::npos) {
			if (pos0 == 0) {
				s.clear();
			} else {
				s = s.substr(0, pos0 - 1);
			}
		}
		if ((pos0 = s.find("/*")) != std::string::npos) {
			if (pos0 == 0) {
				s.clear();
			} else {
				s = s.substr(0, pos0 - 1);
			}
		}
		if ((pos0 = s.find("\\")) != std::string::npos) {
			if (pos0 == 0) {
				s.clear();
			} else {
				s = s.substr(0, pos0 - 1);
			}
		}
		//Removed to integrate Constraint Names
		//if ((pos0 = s.find(":")) != std::string::npos) {
		//	s = s.substr(pos0+1);
		//}
		return s;
	}

	static double round(double r, int tail) {
		if (r > 0.0)
			return floor(r * pow(10.0, tail) + 0.5) * pow(10.0, -tail);
		return ceil(r * pow(10.0, tail) - 0.5) * pow(10.0, -tail);
	}

	static void PAUSE() {
		char t;
		std::cout << "Press a key then press enter (q to quit): ";
		std::cin >> t;
		if (t == 'q')
			exit(0);
	}

	static bool PRESS_KEY(char c) {
		char t;
		std::cout << "Enter key: ";
		std::cin >> t;
		return (c == t);
	}

	static void printLines(const std::list<std::string>& list) {
		std::list<std::string>::const_iterator it = list.begin();
		std::list<std::string>::const_iterator end = list.end();
		for (; it != end; ++it) {
			std::cout << *it << std::endl;
		}
	}
};
}
#endif /*TOOLBOX_HPP_*/
