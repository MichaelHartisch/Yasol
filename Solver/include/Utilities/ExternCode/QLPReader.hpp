/*
*
* Solver: QLPReader.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPREADER_HPP
#define QLPREADER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <list>
#include <vector>
#include <stdexcept>
//#include <gmpxx.h>

#include "Datastructures/Datastructures.hpp"



class QLPReader
{

	public:
		static data::Qlp read( const char* filename, const std::string randPrefix = "", const std::string allPrefix = "" );

	private:

		struct readerHelper {
			mpq_class mult;
			int ind;
		};

		struct constraintHelper {
			std::vector<readerHelper> constraint;
			int type;
			mpq_class rhs;
		};

		class stringHelper{

			public:
				stringHelper( std::string str ) :
					str( str ),
					begin( str.c_str() ),
					end( begin + str.length() * sizeof( char ) ){}

				bool removeFromBeg( unsigned int n ){
					if( begin >= end ){
						return false;
					}

					begin += n * sizeof( char );

					return true;
				}

				bool removeFromEnd( unsigned int n ){
					if( begin >= end ){
						return false;
					}

					end -= n * sizeof( char );

					return true;
				}

				void trimLeft(){
					while( begin < end && ( *begin == ' ' || *begin == '\n' || *begin == '\t' || *begin == '\r' ) ){
						++begin;
					}
				}

				void trimRight(){
					while( begin < end && ( *end == ' ' || *end == '\n' || *end == '\t' || *end == '\r' ) ){
						--end;
					}
				}

				void trim(){
					trimLeft();
					trimRight();
				}

				bool removeUntilFirstOrDoNothing( char c ){
					const char* foo = begin;
					while( foo < end ){
						if( *foo == c ){
							begin = foo + sizeof( char );
							return true;
						}
						++foo;
					}
					return false;
				}

				unsigned int length(){
					return ( end - begin ) / sizeof( char );
				}

				const char& operator[]( unsigned int i ){
					if( begin + sizeof( char ) * i >= end ){
						throw std::runtime_error( "Index out of bounds." );
					}
					return *( begin + sizeof( char ) * i );
				}

				std::string substringUntilFirstOfAndShrink( std::string needles ){
					if( begin >= end ){
						return std::string();
					}
					const char* foo = begin;
					for( ; foo < end; ++foo ){
						bool stop = false;
						for( unsigned int i = 0; i < needles.size(); ++i ){
							if( needles[i] == *foo ){
								++foo;
								stop = true;
								break;
							}
						}
						if( stop ){
							break;
						}
					}

					std::string ret = std::string( begin, ( foo - begin ) / sizeof( char ) );

					begin = foo;
					return ret;
				}


				std::string value(){
					return std::string( begin, ( end - begin ) / sizeof( char ) );
				}

			private:
				const std::string str;
				const char* begin;
				const char* end;
		};

		QLPReader();
		static void parseConstraint( const std::string& s, const std::map<std::string, int>& varNumbers, std::vector<constraintHelper>& matrix );
		static void readVariableAndShrink( stringHelper& s, const std::map<std::string, int>& varNumbers, int& n );
		static std::string trim( const std::string& s );
		static std::string removeComments( const std::string& s );
		static std::string strtolower( const std::string& s );
		static const bool isNumeric( const std::string& s );
		static void readNumberAndShrink( stringHelper& s, mpq_class& n );



//		template<typename T> static void showVector( const std::vector<T> v, const std::string name = "" ) {
//			if (!name.empty()) {
//				std::cout << name << ": ";
//			}
//			std::cout << "{ ";
//			for (unsigned int i = 0; i < v.size(); ++i) {
//				std::cout << v[i] << " ";
//			}
//			std::cout << "}" << std::endl;
//		}

};



#endif
