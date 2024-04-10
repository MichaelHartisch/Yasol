/*
 *
 * Solver: Parser.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Settings/Settings.hpp"
#include "Utilities/Parser.hpp"
#include "Utilities/ToolBox.hpp"
#include "Utilities/QlpConverter.hpp"
//#include "QlpConverter.hpp"
using namespace std;
namespace utils {
std::string Parser::LOG_TAG = "Parser";

void Parser::createQlp(const std::string& path, data::Qlp& p) {
    // Clear the target
    p.clear();
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Reading File... : " + path);
    // 1. Create a Parser Instance and read the input file
    Parser parser(path);
    // 2. Sort all lines by their membership and add them to special string lists
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Sorting Lines...");
    parser.sortLines();
    // 3. Create all columns of the Qlp and the objective function
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Columns...");
    parser.createQlpColumns();
    // 4. Create all rows of the Qlp
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Rows...");
    parser.createRows();
    // 5. Create qlp
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating QLP...");
    parser.createQlp(p);
    // 6. QLPs are bounded above an below
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Setting Bounds...");
    std::vector<data::QpVar*> v = p.getVariableVector();
    for (unsigned int i = 0; i < v.size(); i++) {
        if (v[i]->getLowerBound() < LB_QLP)
            v[i]->setLowerBound(LB_QLP);
        if (v[i]->getUpperBound() > UB_QLP)
            v[i]->setUpperBound(UB_QLP);
    }
    
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Done");
    std::vector<data::QpRhs> rhs = p.getRhsVec();
    //std::cerr << "in create: ";
    //for (int i = 0; i < rhs.size();i++) {
    //    std::cerr << rhs[i].getResponsibility();
    //}
    //std::cerr << std::endl;
}

void Parser::createQlp(const std::string& path, data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs) {
    //Clear target
    obj.clear();
    vars.clear();
    matrix.clear();
    rhs.clear();
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Reading File... : " + path);
    // 1. Create a Parser Instance and read the input file
    Parser parser(path);
    
    // 2. Sort all lines by their membership and add them to special string lists
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Sorting Lines...");
    parser.sortLines();
    
    // 3. Create all columns of the Qlp and the objective function
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Columns...");
    parser.createQlpColumns(vars, obj);
    
    // 4. Create all rows of the Qlp
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Rows...");
    parser.createRows(vars, matrix, rhs);
    
    // 6. QLPs are bounded above an below
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Setting Bounds...");
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (vars[i].getLowerBound() < LB_QLP)
            vars[i].setLowerBound(LB_QLP);
        if (vars[i].getUpperBound() > UB_QLP)
            vars[i].setUpperBound(UB_QLP);
    }
    
    // 7. Make the problem an minimization problem
    if (obj.getObjective() == data::QpObjFunc::max)
        obj.reverseObjectiveFuntion();
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Done");
    parser.clear();
    //std::cerr << "in create 2: ";
    //for (int i = 0; i < rhs.size();i++) {
    //    std::cerr << rhs[i].getResponsibility();
    //}
    //std::cerr << std::endl;
    
}

void Parser::createQlp(data::Qlp& qlp) {
    qlp.initWithQlpParts(_obj,_vars,_lhsIE,_rhs);
    unsigned int size = qlp.getConstraintCount();
    //utils::QlpConverter::removeEmptyConstraints(qlp);
    if(size!=qlp.getConstraintCount()){
        if (LOG_PARSER)
            utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "EMPTY CONSTRAINT REMOVED...");
    }
    //std::cerr << "in create 3: ";
    //std::vector<data::QpRhs> rhs = qlp.getRhsVec();
    //for (int i = 0; i < rhs.size();i++) {
    //    std::cerr << rhs[i].getResponsibility();
    //}
    //std::cerr << std::endl;
    //std::vector<data::QpRhs> rhs = qlp.getRhsVec();
    //for (int i = 0; i < _rhs.size();i++) {
    //    std::cerr << _rhs[i].getResponsibility();
    //}
    //std::cerr << std::endl;
}

Parser::Parser(const std::string& inFile) :
objective(data::QpObjFunc::min), qlp(), objectiveFunction(), constraints(), u_constraints(),  bounds(), generals(), binaries(), all(), exists(), random(), order(), _obj(), _vars(), _rhs(), _lhsIE(), nameToVariableIndex() {
    utils::ToolBox::readFromFile(inFile, qlp);
    std::list<std::string>::iterator start = qlp.begin();
    std::list<std::string>::iterator end = qlp.end();
    while (start != end) {
        utils::ToolBox::removeComments(*start);
        utils::ToolBox::removeStartEndWhitespaces(*start);
        start++;
    }
}

int Parser::getVariableIndexByName(const std::string& varName) const {
    std::unordered_map<std::string, unsigned int>::const_iterator iter;
    if ((iter = nameToVariableIndex.find(varName)) != nameToVariableIndex.end()) {
        return iter->second;
    } else {
        std::cout<<"Name not found: "<<varName<<std::endl;
        throw utils::DataStructureException("Parser: VarName not found...");
        return -1;
    }
}

data::QpVar& Parser::getVarByName(const std::string& name, std::vector<
                                  data::QpVar>& vec) const {
    int index = this->getVariableIndexByName(name);
    if(index>=0&&index<vec.size()){
        return vec[index];
    }else{
        utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Name: "
                                 + name);
        utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable",
                                 "Variables: " + data::QpVar::vecToString(vec));
        throw utils::DataStructureException(
                                            "getQpVarByName --> variable does not exist: "+name);
    }}


void Parser::printLines(std::list<std::string>& list) {
    std::list<std::string>::iterator it = list.begin();
    std::list<std::string>::iterator end = list.end();
    for (; it != end; ++it) {
        cout << *it << endl;
    }
}

void Parser::createQlpColumns() {
    this->createQlpColumns(this->_vars, this->_obj);
    if (this->_vars.back().getQuantifier() != data::QpVar::exists){
        this->_vars.push_back(data::QpVar("AuxVar", this->_vars.size(), 0, 0));
        this->nameToVariableIndex.insert(std::pair<std::string, int>("AuxVar", _vars.size()-1));
        this->getVarByName("AuxVar", _vars).setQuantifier(data::QpVar::exists);
    }
}
void Parser::createRows() {
    this->createRows(this->_vars, this->_lhsIE, this->_rhs);
}

void Parser::createQlpColumns(std::vector<data::QpVar>& v, data::QpObjFunc& o) {
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating variable order...");
    createVariableOrder(v);
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating quantifiers...");
    addQuantifiers(v);
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating integrals...");
    addIntegrals(v);
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating bounds...");
    addBounds(v);
    if (LOG_PARSER)
        utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating objective function...");
    addObjectiveFunction(v, o);
}

void Parser::createRows(std::vector<data::QpVar>& vars, data::QpSparseMatrix& lhsIE, std::vector<data::QpRhs>& rhs) {
    createRows(vars,lhsIE,rhs,"EXISTS");
    createRows(vars,lhsIE,rhs,"UNIVERSAL");
}

void Parser::CheckDigits(std::string& Coefficient, bool IsScientific){
    int MoreSymbols=0;
    string::size_type position= Coefficient.find(".");
    if(position != std::string::npos){
        MoreSymbols++;
        if(IsScientific && position != 1) throw utils::ParserException("Have the separator at second position in the coefficient when using scientific notation " + Coefficient);
        else if (IsScientific && position == 1 && Coefficient[0]=='0')  throw utils::ParserException("Use a number between 1 and 10 for the coefficient when using scientific notation " + Coefficient);
    }
    else if(IsScientific && (Coefficient.size()>1 || Coefficient[0]=='0')) throw utils::ParserException("Use a number between 1 and 10 for the coefficient when using scientific notation " + Coefficient);
    if(Coefficient.find("-") != std::string::npos || Coefficient.find("+") != std::string::npos) MoreSymbols++;
    if(Coefficient.size()>MAX_DIGITS+MoreSymbols)
        throw utils::ParserException("Do not use coefficients with more than " +std::to_string(MAX_DIGITS) +" digits: " + Coefficient);
}

bool Parser::BoundIsInf(std::string& bound){
    std::string valueStr, exponentStr;
    char c;
    double sign = 1;
    data::QpNum AuxSciCoef; //Variable for quickly checking scientific Notation
    data::QpNum AuxSciExp; //Variable for quickly checking scientific Notation
    bool scientificMode = false;
    bool FoundFirstSign = false;
    valueStr="";
    exponentStr="";
    
    for (unsigned int i = 0; i < bound.size(); i++) {
        if((c = bound[i]) != '.' && !isdigit(c) && c != ' '){
            if (c=='+' || c=='-'){
                if(FoundFirstSign) // If the +/- is not the first one found (scientifc +/- is scipped right away as soon as e/E is found)
                    throw utils::ParserException("only numerical values allowed for bound of variabe; no calculations permitted " + bound);
                if (!FoundFirstSign){
                    if(c=='-') sign = -1;
                }
            }
            else if ((c=='e' || c=='E') &&      //Found potential scientific mode
                     !scientificMode &&              //And we are not yet in scientifc mode
                     i+1 < bound.size() && //There is still something right of e/E
                     (bound[i+1]=='-' || bound[i+1]=='+' || isdigit(bound[i+1]))) //which is some number
            {
                CheckDigits(valueStr,1);
                AuxSciCoef=valueStr;
                scientificMode=true;
                valueStr+=c;
                if(bound[i+1]=='-' || bound[i+1]=='+'){
                    exponentStr+=bound[i+1];
                    valueStr+=bound[i+1];
                    i++;    //skip explicitly reading the sign; but add it to the strings to be handeled by QpRHS/QpVal
                    continue;
                }
            }
            else throw utils::ParserException("non-numerical at bounds " + bound);
        }
        if (c=='.' && scientificMode) throw utils::ParserException("only use integers as exponents for scientific notation " + bound);
        else if (c=='.') valueStr+=c;
        if (isdigit(c)){
            if (scientificMode)
                exponentStr+=c; //exponentStr for quickly checking the size
            valueStr+=c;
        }
        if (c!=' ' && !FoundFirstSign) FoundFirstSign=true;
    }
    
    bool IsInf=false;
    string::size_type position= valueStr.find(".");
    if (scientificMode){
        if(position != std::string::npos && position != 1)
            throw utils::ParserException("Have the separator at second position when using scientific notation " + valueStr);
        AuxSciExp=exponentStr;
        IsInf=!CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_BOUND);
        if(AuxSciCoef==0 && AuxSciExp==0) bound="0";
    }
    else{
        CheckDigits(valueStr);
        if(position != std::string::npos){
            if(valueStr.substr(0,position).size() >MAX_EXP_BND) IsInf=true;
        }
        
        else{
            if(valueStr.size()>MAX_EXP_BND) IsInf=true;
        }
        //double val=std::stod(valueStr);
        //if (abs(val)>=pow(10,15))
        //   IsInf=true;
    }
    if(IsInf) std::cerr << "Warning: Variable is interpreted as unbounded. Bound is set to " << UB_QLP << " in case of upper bound and " << LB_QLP << " in case of a lower bound." << std::endl;
    return IsInf;
    
}
bool Parser::CheckScientificNotation(data::QpNum& coefficient, data::QpNum& exponent, const int Mode){
    switch (Mode) {
        case ST_RHS:
            if (exponent >= MAX_EXP_RHS){
                throw utils::ParserException("Exponent of scientific RHS is too large: " + std::to_string(coefficient.asDouble()) + "e"+std::to_string(exponent.asDouble()));
            }
            else if (exponent <= MIN_EXP_RHS){
                std::cerr << "Warning: exponent of scientific rhs is too small: " << std::to_string(coefficient.asDouble()) << "e"<<std::to_string(exponent.asDouble())<<". Rhs is set to zero."<<std::endl;
                coefficient = 0;
                exponent = 0;
            }
            return true;
        case ST_BOUND:
            if (exponent >= MAX_EXP_BND){
                std::cerr << "Warning: Variable Bound is too large. Scientific notation was used with exponent>=15. Bound is reduced to "<< MAX_EXP_BND << " and we proceed.";
                coefficient = 1;
                exponent = 15;
            }
            else if (exponent <= MIN_EXP_BND){
                std::cerr << "Warning: Exponent of scientific bound is too small: " << std::to_string(coefficient.asDouble()) << "e"<<std::to_string(exponent.asDouble())<< ". Bound is set to zero."<<std::endl;
                coefficient = 0;
                exponent = 0;
            }
            return true;
            break;
        case ST_COEF:
            if (exponent >= MAX_EXP_COEF){
                throw utils::ParserException("Exponent of scientific coefficient is too large: " + std::to_string(coefficient.asDouble()) + "e"+std::to_string(exponent.asDouble()));
            }
            else if (exponent <= MIN_EXP_COEF){
                std::cerr << "Warning: Exponent of scientific coefficient is too small: " << std::to_string(coefficient.asDouble()) << "e"<<std::to_string(exponent.asDouble()) << ". Coefficient is set to zero."<< std::endl;
                coefficient = 0;
                exponent = 0;
            }
            return true;
        default:
            throw "Unsupported Mode for Checking";
            return false;
    }
}
void Parser::createRows(std::vector<data::QpVar>& vars, data::QpSparseMatrix& lhsIE, std::vector<data::QpRhs>& rhs,std::string Q) {
    std::string ratioSign, valueStr, varName, constraint, points, constraintName, exponentStr;
    std::string::size_type separator;
    bool DetectedDeprecatedNamingScheme=false;
    std::list<std::string>::const_iterator it;
    std::list<std::string>::const_iterator end;
    if(Q.compare("EXISTS")!=0){
        it = constraints.begin();
        end = constraints.end();
    }
    else if(Q.compare("UNIVERSAL")!=0){
        it = u_constraints.begin();
        end = u_constraints.end();
    }
    else throw utils::ParserException("Unknown Quantifier " + Q);

    char c;
    double sign = 1;
    data::QpNum val;
    data::QpNum AuxSciCoef; //Variable for quickly checking scientific Notation
    data::QpNum AuxSciExp; //Variable for quickly checking scientific Notation
    data::QpNum rhsVal;
    bool valueMode = false;
    bool varNameMode = false;
    bool scientificMode = false;
    bool FoundFirstSignRHS = false;
    unsigned int index;
    
    for (int i = 0; it != end; it++, ++i) {
        
        constraintName.clear();
        
        constraint = *it;
        
        if ((separator = constraint.find(":")) != std::string::npos) {
            constraintName = constraint.substr(0,separator);
            //cerr << "CONSTRAINTNAME:" << constraintName << endl;
            constraint = constraint.substr(separator+1);
            utils::ToolBox::removeStartEndWhitespaces(constraintName);
            utils::ToolBox::removeStartEndWhitespaces(constraint);
            //std::cout << constraintName << std::endl;
        }
        
        if ((separator = constraint.find("=")) == string::npos) {
            throw utils::ParserException("Unknown RatioSign at constraint " + constraint);
        }
        else{
            //NEW ATTENTION: Explicitly check RHS side for Variables, i.e. non-number entries, and throw an error
            //Also new: Collect RHS and check for scientific notation
            valueStr = "";
            sign = 1;
            exponentStr = "";
            
            scientificMode=false;
            FoundFirstSignRHS=false;
            
            for (unsigned int i = 0; i < constraint.substr(separator + 1).size(); i++) {
                if((c = constraint.substr(separator + 1)[i]) != '.' && !isdigit(c) && c != ' '){
                    if (c=='+' || c=='-'){
                        if(FoundFirstSignRHS) // If the +/- is not the first one found (scientifc +/- is scipped right away as soon as e/E is found)
                            throw utils::ParserException("only numerical values allowed at rhs of constraints; no calculations permitted " + constraint);
                        if (!FoundFirstSignRHS){
                            if(c=='-') sign = -1;
                        }
                    }
                    else if ((c=='e' || c=='E') &&      //Found potential scientific mode
                             !scientificMode &&              //And we are not yet in scientifc mode
                             i+1 < constraint.substr(separator + 1).size() && //There is still something right of e/E
                             (constraint.substr(separator + 1)[i+1]=='-' || constraint.substr(separator + 1)[i+1]=='+' || isdigit(constraint.substr(separator + 1)[i+1]))) //which is some number
                    {
                        CheckDigits(valueStr,1);
                        AuxSciCoef=valueStr;
                        scientificMode=true;
                        valueStr+=c;
                        if(constraint.substr(separator + 1)[i+1]=='-' || constraint.substr(separator + 1)[i+1]=='+'){
                            exponentStr+=constraint.substr(separator + 1)[i+1];
                            valueStr+=constraint.substr(separator + 1)[i+1];
                            i++;    //skip explicitly reading the sign; but add it to the strings to be handeled by QpRHS/QpVal
                            continue;
                        }
                    }
                    else throw utils::ParserException("non-numerical at rhs of constraints " + constraint);
                }
                else if (c == '.'){
                    if (scientificMode) throw utils::ParserException("only use integers as exponents for scientific notation " + constraint);
                    valueStr+=c;
                }
                if (isdigit(c)){
                    if (scientificMode)
                        exponentStr+=c; //exponentStr for quickly checking the size
                    valueStr+=c;
                }
                if (c!=' ' && !FoundFirstSignRHS) FoundFirstSignRHS=true;
            }
            
            if (scientificMode){
                AuxSciExp=exponentStr;
                CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_RHS);
                if(AuxSciCoef==0 && AuxSciExp==0) valueStr="0";
                AuxSciCoef=0;
                AuxSciExp=0;
                if(exponentStr.size()==1 && (exponentStr[0]=='+' || exponentStr[0]=='-'))
                    throw utils::ParserException("Unsupported scientific notation. Did you use e/E as variable name? " + constraint);
                exponentStr="";
                scientificMode = false;
            }
            else CheckDigits(valueStr);
            rhsVal=valueStr;
            rhsVal*=sign;
            if ((separator = constraint.find("<=")) != string::npos) {
                rhs.push_back(data::QpRhs(rhsVal,data::QpRhs::smallerThanOrEqual));
            } else if ((separator = constraint.find(">=")) != string::npos) {
                rhs.push_back(data::QpRhs(rhsVal,data::QpRhs::greaterThanOrEqual));
            } else if ((separator = constraint.find("=")) != string::npos) {
                rhs.push_back(data::QpRhs(rhsVal,data::QpRhs::equal));
            }
        }
        
        
        if(Q.compare("EXISTS")!=0){
            if (constraintName.size() > 1 && constraintName[0] == 'U' && constraintName[1] == '_') {
                if(!DetectedDeprecatedNamingScheme){
                    cerr << "WARNING: Detected deprecated naming scheme for universal constraint." << endl;
                    cerr << "Constraints within SUBJECT TO, with names starting with \"U_\" are (still) considered universal constraints." << endl;
                    cerr << "Better use the keyword UNCERTAINTY SUBJECT TO instead." << endl;
                    DetectedDeprecatedNamingScheme=true;
                }
                rhs[rhs.size()-1].setResponsibility(data::Constraint::UNIVERSAL);
            } else {
                //cerr << "detected THISORTHAT constrant at Index" << rhs.size()-1 << endl;
                rhs[rhs.size()-1].setResponsibility(data::Constraint::EXISTENTIAL);
            }
        }
        else{
            rhs[rhs.size()-1].setResponsibility(data::Constraint::UNIVERSAL);
        }
        //std::cerr << "Control:";
        //for (int j = 0; j < rhs.size();j++) {
        //    std::cerr << rhs[j].getResponsibility();
        //}
        //cerr << ";";
        
        lhsIE.push_back(std::vector<data::IndexedElement>());
        index = lhsIE.size() - 1;
        valueStr = "1";
        exponentStr = "";
        sign = 1;
        points = constraint.substr(0, separator);
        
        for (unsigned int i = 0; i < points.size(); i++) {
            // skip white spaces
            if ((c = points[i]) == ' ') {
            } else if (c == '+' || c == '-') {
                // we read a plus or minus sign, hence we reached a new point
                // NEW: UNLESS there is a number in scientific notation
                if (valueMode) {
                    // error because there is no variable, only a value
                    throw utils::ParserException("Unsupported Point mode");
                } else if (varNameMode) {
                    
                    if (scientificMode){
                        AuxSciExp=exponentStr;
                        CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_COEF);
                        if(AuxSciCoef==0 && AuxSciExp==0) valueStr="0";
                        AuxSciCoef=0;
                        AuxSciExp=0;
                        if(exponentStr.size()==1 && (exponentStr[0]=='+' || exponentStr[0]=='-'))
                            throw utils::ParserException("Unsupported scientific notation. Did you use e/E as variable name? " + constraint);
                        exponentStr="";
                        scientificMode = false;
                    }
                    else CheckDigits(valueStr);
                    val = valueStr;
                    
                    lhsIE[index].push_back(data::IndexedElement(getVariableIndexByName(varName), data::QpNum(val *sign)));
                    //std::cerr << varName << " " << val.asDouble() << " " << sign << " " << shift.asDouble() << endl;
                    
                    valueMode = false;
                    varNameMode = false;
                }
                // the new sign
                c == '+' ? sign = 1 : sign = -1;
                valueStr = "1";
            } else if (isdigit(c) || c == '.') {
                if (valueMode) {
                    valueStr += c;
                    if (scientificMode){
                        if(c == '.')
                            throw utils::ParserException("only use integers as exponents for scientific notation " + constraint);
                        exponentStr += c;
                    }
                }
                else if (varNameMode) {
                    varName += c;
                } else {
                    valueMode = true;
                    valueStr = c;
                }
            } else /*if (isalnum(c) || c == '_')*/{
                if (varNameMode) {
                    varName += c;
                } else if (valueMode || scientificMode) {
                    // end value mode and start verNameMode
                    if((c == 'e' || c == 'E')
                        && i+1 < points.size() && (points[i+1] == '-' || points[i+1] =='+' || isdigit(points[i+1]))){
	                bool PossibleVariable=false;
			if(isdigit(points[i+1])){
			    //Maybe this is in fact a variable name like E1 or e2411; Note that variable names like e1d are forbidden
			    PossibleVariable=true;
			    int loc=i+1;
			    while(loc < points.size() && points[loc]!= '+' && points[loc]!= '-' && points[loc]!= ' ' && points[loc]!= '<' && points[loc]!= '>' && points[loc]!= '='){
			        if(isdigit(points[loc])) loc++;
			        else{
			     	    PossibleVariable=false;
				    break;
				}
		            }
			}
			if(!PossibleVariable){
                            if (scientificMode)
                                throw utils::ParserException("Unsupported Double Exponent in scientific notation. Also do not use e/E as variable name!");
                            CheckDigits(valueStr,1);
                            AuxSciCoef=valueStr;
                            scientificMode=true;
                            valueStr+=c;
                            if(points[i+1]=='-' || points[i+1]=='+'){
                                exponentStr+=points[i+1];
                                valueStr+=points[i+1];
                                i++;    //skip explicitly reading the sign of exponent; but add it to the strings to be handeled by QpRHS/QpVal
                                continue;
                            }
                        }
                        else{//Found possible variable like e121 or E7
			    varName = c;
                            valueMode = false;
                            varNameMode = true;
                        }
                    }
                    else{
                        varName = c;
                        valueMode = false;
                        varNameMode = true;
                    }
                } else {
                    varName = c;
                    varNameMode = true;
                }
            }
            if (i == points.size() - 1) {
                // we reached the end of the point string for this row
                /*val = valueStr;
                 lhsIE[index].push_back(data::IndexedElement(getVariableIndexByName(varName), data::QpNum(val * sign)));
                 sign = 1;
                 valueStr = "1";
                 valueMode = false;
                 varNameMode = false;*/
                
                if (scientificMode){
                    AuxSciExp=exponentStr;
                    CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_COEF);
                    if(AuxSciCoef==0 && AuxSciExp==0) valueStr="0";
                    AuxSciCoef=0;
                    AuxSciExp=0;
                    if(exponentStr.size()==1 && (exponentStr[0]=='+' || exponentStr[0]=='-'))
                        throw utils::ParserException("Unsupported scientific notation. Did you use e/E as variable name? " + constraint);
                    exponentStr="";
                    scientificMode = false;
                }
                else CheckDigits(valueStr);
                val = valueStr;
                
                lhsIE[index].push_back(data::IndexedElement(getVariableIndexByName(varName), data::QpNum(val *sign)));
                sign = 1;
                valueStr = "1";
                valueMode = false;
                varNameMode = false;
            }
        }
    }
    //std::cerr << "in create 0: ";
    //for (int i = 0; i < rhs.size();i++) {
    //    std::cerr << rhs[i].getResponsibility();
    //}
    //std::cerr << std::endl;
}

void Parser::sortLines() {
    
    int mode = Parser::INIT;
    string line, currentConstraint;
    string::size_type pos0;
    list<string>::const_iterator it = qlp.begin();
    list<string>::const_iterator en = qlp.end();
    
    std::string end("end"),End("End"),END("END");
    
    for (; it != en; ++it) {
        if ((line = *it).size() > 0 && line.compare(end) && line.compare(End) && line.compare(END)) {
            // check which mode is the current active one
            switch (mode) {
                case INIT:
                    // determine the objective
                    if (line == "MINIMIZE" || line == "Minimize" || line == "Min" || line == "MIN") {
                        objective = data::QpObjFunc::min;
                    } else if (line == "MAXIMIZE" || line == "Maximize" || line == "Max" || line == "MAX") {
                        objective = data::QpObjFunc::max;
                    } else
                        throw utils::ParserException("Unrecognized keyword in line: " + line);
                    // change the mode since we read the max or min keyword
                    mode = OBJECTIVE_FUNCTION;
                    break;
                case OBJECTIVE_FUNCTION:
                    if (line == "SUBJECT TO" || line == "Subject to" ||line == "subject to" || line == "Subject To" || line == "s.t." || line == "st." || line == "st") {
                        // change the mode
                        mode = CONSTRAINTS;
                    } else {
                        // remove the description of the objective function and trim the result
                        if ((pos0 = line.find(":")) != string::npos) {
                            utils::ToolBox::removeStartEndWhitespaces(line = line.substr(pos0 + 1));
                        }
                        // add the line to the string that contains the whole objective function
                        objectiveFunction += line;
                    }
                    break;
                case CONSTRAINTS:
                    if (line == "BOUNDS" || line == "Bounds" || line == "bounds") {
                        // change the mode
                        mode = BOUNDS;
                    } else if (line == "UNCERTAINTY SUBJECT TO" || line == "u s.t." || line == "u.s.t." || line == "ust"  ||  line == "UNCERTAINTY Subject to" || line == "Uncertainty Subject to" || line == "uncertainty Subject to" ||line == "UNCERTAINTY subject to"||line == "Uncertainty subject to" ||line == "uncertainty subject to"  || line == "UNCERTAINTY Subject To" || line == "Uncertainty Subject To" ||line == "uncertainty Subject To" ||line == "UNCERTAINTY s.t." ||  line == "Uncertainty s.t." ||  line == "uncertainty s.t." || line == "UNCERTAINTY st." || line == "Uncertainty st." || line == "uncertainty st." ||  line == "uncertainty st" || line == "UNCERTAINTY st" || line == "Uncertainty st") {
                        mode = UNCERTAINTY;
                    }else if (line == "EXISTS" || line == "EXISTS" || line == "Exists" || line == "exists" || line == "Exist" || line == "exist") {
                        mode = EXISTS;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        if ((pos0 = line.find(":")) != string::npos) {
                            // this is a new constraint with a comment
                            //if(utils::ToolBox::removeStartEndWhitespaces(line = line.substr(pos0 + 1)).size());
                            //	currentConstraint = line;
                        }
                        if ((pos0 = line.find("=")) != string::npos) {
                            // finish the previous constraint and add it to the list
                            currentConstraint+=line;
                            constraints.push_back(currentConstraint);
                            currentConstraint.clear();
                        } else {
                            currentConstraint += line;
                        }
                        
                    }
                    break;
                case UNCERTAINTY:
                    if (line == "BOUNDS" || line == "Bounds" || line == "bounds") {
                        // change the mode
                        mode = BOUNDS;
                    } else if (line == "EXISTS" || line == "EXISTS" || line == "Exists" || line == "exists" || line == "Exist" || line == "exist") {
                        mode = EXISTS;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        if ((pos0 = line.find(":")) != string::npos) {
                            //line="U_"+line;
                            // this is a new constraint with a comment
                            //if(utils::ToolBox::removeStartEndWhitespaces(line = line.substr(pos0 + 1)).size());
                            //      currentConstraint = line;
                        }
                        //else line="U_1: "+line;
                        if ((pos0 = line.find("=")) != string::npos) {
                            // finish the previous constraint and add it to the list
                            currentConstraint+=line;
                            u_constraints.push_back(currentConstraint);
                            currentConstraint.clear();
                        } else {
                            currentConstraint += line;
                        }
                    }
                    break;
                case BOUNDS:
                    if (line == "GENERALS" || line == "Generals" || line == "generals"|| line == "GENERAL" || line == "General" || line == "general" || line == "GEN" || line == "GENERAL") {
                        mode = GENERALS;
                    } else if (line == "BINARIES" || line == "Binaries" || line == "binaries" || line == "BIN" || line == "Bin" || line == "bin" || line == "BINARY" || line == "Binary" || line == "binary") {
                        mode = BINARIES;
                    } else if (line == "EXISTS" || line == "Exists" || line == "exists" || line == "Exist" || line == "exist") {
                        mode = EXISTS;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        bounds.push_back(line);
                    }
                    break;
                case GENERALS:
                    if (line == "BINARIES" || line == "Binaries" || line == "BIN" || line == "BINARY") {
                        mode = BINARIES;
                    } else if (line == "EXISTS" || line == "Exists" || line == "exists" || line == "Exist" || line == "exist") {
                        mode = EXISTS;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        generals += " " + line;
                    }
                    break;
                case BINARIES:
                    if (line == "GENERALS" || line == "Generals" || line == "generals"|| line == "GENERAL" || line == "General" || line == "general" || line == "GEN" || line == "GENERAL") {
                        mode = GENERALS;
                    } else if (line == "EXISTS" || line == "Exists" || line == "exists" || line == "Exist" || line == "exist") {
                        mode = EXISTS;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        binaries += " " + line;
                    }
                    break;
                case EXISTS:
                    if (line == "ORDER" || line == "Order") {
                        mode = ORDER;
                    } else if (line == "ALL" || line == "All" || line == "all") {
                        mode = ALL;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        exists += " " + line;
                    }
                    break;
                case ALL:
                    if (line == "ORDER" || line == "Order") {
                        mode = ORDER;
                    } else if (line == "RANDOM" || line == "Random") {
                        mode = RANDOM;
                    } else {
                        all += " " + line;
                    }
                    break;
                case RANDOM:
                    if (line == "ORDER" || line == "Order") {
                        mode = ORDER;
                    } else {
                        random += " " + line;
                    }
                    break;
                case ORDER:
                    order += " " + line;
                    break;
                default:
                    throw "Unsupported Quantifier";
            }
        }
    }
    qlp.clear();
}

void Parser::addQuantifiers(std::vector<data::QpVar>& vars) {
    string l;
    bool readVar = false;
    // check if existential quantifiers were added to the qlp
    if (exists.size() > 0) {
        for (unsigned int i = 0; i < exists.size(); ++i) {
            // skip white spaces
            if (exists[i] == ' ' && !readVar) {
            } else if ((exists[i] == ' ' && readVar) || i == exists.size() - 1) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == exists.size() - 1)
                    l += exists[i];
                this->getVarByName(l, vars).setQuantifier(data::QpVar::exists);
                l = "";
                readVar = false;
            } else {
                readVar = true;
                // the current char is part of the name of the variable
                l += exists[i];
            }
        }
    }
    
    if (all.size() > 0) {
        for (unsigned int i = 0; i < all.size(); ++i) {
            if (all[i] == ' ' && !readVar) {
            } else if ((all[i] == ' ' && readVar) || i == all.size() - 1) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == all.size() - 1)
                    l += all[i];
                this->getVarByName(l, vars).setQuantifier(data::QpVar::all);
                l = "";
                readVar = false;
            } else {
                readVar = true;
                // the current char is part of the name of the variable
                l += all[i];
            }
        }
    }
    
    if (random.size() > 0) {
        for (unsigned int i = 0; i < random.size(); ++i) {
            if ((random[i] == ' ') && !readVar) {
            } else if ((random[i] == ' ' && readVar) || i == random.size() - 1) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == random.size() - 1)
                    l += random[i];
                this->getVarByName(l, vars).setQuantifier(data::QpVar::random);
                l = "";
                readVar = false;
            } else {
                readVar = true;
                // the current char is part of the name of the variable
                l += random[i];
            }
        }
    }
}

void Parser::addIntegrals(std::vector<data::QpVar>& vars) {
    string l;
    bool readVar = false;
    if (generals.size() > 0) {
        for (unsigned int i = 0; i < generals.size(); ++i) {
            // skip white spaces
            if (generals[i] == ' ' && !readVar) {
            } else if ((generals[i] == ' ' && readVar) || i == generals.size() - 1) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == generals.size() - 1)
                    l += generals[i];
                this->getVarByName(l, vars).setNumberType(data::QpVar::generals);
                l = "";
                readVar = false;
            } else {
                readVar = true;
                // the current char is part of the name of the variable
                l += generals[i];
            }
        }
    }
    if (binaries.size() > 0) {
        for (unsigned int i = 0; i < binaries.size(); ++i) {
            // skip white spaces
            if (binaries[i] == ' ' && !readVar) {
            } else if ((binaries[i] == ' ' && readVar) || i == binaries.size() - 1) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == binaries.size() - 1)
                    l += binaries[i];
                data::QpVar& v = this->getVarByName(l, vars);
                v.setNumberType(data::QpVar::binaries);
                v.setBounds(0, 1);
                l = "";
                readVar = false;
            } else {
                readVar = true;
                // the current char is part of the name of the variable
                l += binaries[i];
            }
        }
    }
}

void Parser::addBounds(std::vector<data::QpVar>& vars) {
    string l, tmp, tmp1, tmp2;
    if (bounds.size() > 0) {
        string bound, lb, varName, ub;
        string::size_type firstSeparator, lastSeparator;
        list<string>::const_iterator it = bounds.begin();
        list<string>::const_iterator end = bounds.end();
        
        for (unsigned int i = 0; it != end; ++it, i++) {
            bound = *it;
            // get positions of first and last separator
            firstSeparator = bound.find("<=");
            lastSeparator = bound.find("<=", firstSeparator + 1);
            // if no smaller than equality sign exists, we look for an equal sign
            if (firstSeparator == string::npos && lastSeparator == string::npos) {
                firstSeparator = lastSeparator = bound.find("=");
            }
            if (firstSeparator < lastSeparator) {
                // lb <= x <= ub
                lb=bound.substr(0, firstSeparator);
                utils::ToolBox::removeStartEndWhitespaces(varName = bound.substr(firstSeparator + 2, lastSeparator - firstSeparator - 2));
                ub = bound.substr(lastSeparator + 2);
            } else if (firstSeparator == lastSeparator && firstSeparator != string::npos) {
                utils::ToolBox::removeStartEndWhitespaces(varName = bound.substr(0, firstSeparator));
                tmp = bound.substr(firstSeparator + 1);
                firstSeparator = tmp.find_first_of("[");
                lastSeparator = tmp.find_last_of("[");
                
                if(firstSeparator == string::npos){
                    lb=ub=utils::ToolBox::removeStartEndWhitespaces(tmp);
                }else{
                    
                    tmp1 = tmp.substr(firstSeparator, lastSeparator);
                    tmp2 = tmp.substr(lastSeparator);
                    std::vector<data::QpNum> vRange;
                    std::vector<data::QpRational> vDistr;
                    data::QpNum::stringVecToQpNumVec(tmp1, vRange);
                    data::QpNum::stringVecToQpNumVec(tmp2, vDistr);
                    this->getVarByName(varName, vars).setVariableRange(vRange, vDistr);
                    tmp.clear();
                    tmp1.clear();
                    tmp2.clear();
                    continue;
                }
                
            } else {
                // something is wrong
                throw utils::ParserException("Invalid Bound: " + bound);
            }
            
            data::QpVar& v = this->getVarByName(varName, vars);
            if(lb == std::string("-inf"))
                v.setUnboundedBelow();
            else if (BoundIsInf(lb))
                v.setUnboundedBelow();
            else
                v.setLowerBound(data::QpNum(lb));
            if(ub == std::string("+inf"))
                v.setUnboundedAbove();
            else if (BoundIsInf(ub))
                v.setUnboundedAbove();
            else
                v.setUpperBound(data::QpNum(ub));
        }
    }
    
}

void Parser::createVariableOrder(std::vector<data::QpVar>& vars) {
    string l;
    bool readVar = false;
    int columnIndex = 0;
    if (order.size() == 0) {
        throw utils::ParserException("No order specified ");
    } else {
        data::QpNum lb(true), ub(false);
        for (unsigned int i = 0, size1 = order.size(), size2 = size1 - 1; i < size1; i++) {
            if (!readVar && order[i] == ' ') {
            } else if ((i == size2) || (readVar && order[i] == ' ')) {
                //we read the current variable to its end and now have to add it to the qlp
                if (i == size2)
                    l += order[i];
		    if((l[0]=='e' || l[0]=='E') &&l.size()==1)
		         throw utils::ParserException("Do not use variable names 'e' or 'E', to avoid ambiguous notation: scientific notation is allowed.");
 		    if((l[0]=='e' || l[0]=='E') && isdigit(l[1])){
                        bool allDigits=true;
                        for (int k=1;k<l.size();k++){
                            if(l[k]==' ') break;
                            else if(!isdigit(l[k])){
                                allDigits=false;
                                break;
                            }
                        }
                        if(!allDigits) throw utils::ParserException("Found variable named " + l +". Do not use variable names starting with 'e' or 'E' preceeded by digits followed by non-digits in order to avoid ambiguous notation: scientific notation is allowed.");
                    }

                vars.push_back(data::QpVar(l, columnIndex, lb, ub));
                this->nameToVariableIndex.insert(std::pair<std::string, int>(l, columnIndex));
                l.clear();
                readVar = false;
                columnIndex++;
            } else {
                readVar = true;
                l += order[i];
            }
        }
    }
}

void Parser::addObjectiveFunction(std::vector<data::QpVar>& vars, data::QpObjFunc& obj) {
    //q.setObjective(objective);
    string ratioSign, valueStr, varName, constraint, points, exponentStr;
    valueStr = "1";
    char c;
    double sign = 1;
    data::QpNum val(1);
    data::QpNum AuxSciCoef; //Variable for quickly checking scientific Notation
    data::QpNum AuxSciExp;
    bool valueMode = false;
    bool varNameMode = false;
    bool scientificMode = false;
    obj.setObjective(objective);
    obj.setSize(vars.size());
    // create the points that correspond to this row
    for (unsigned int i = 0; i < objectiveFunction.size(); i++) {
        // skip white spaces
        if ((c = objectiveFunction[i]) == ' ') {
        } else if (c == '+' || c == '-') {
            // we read a plus or minus sign, hence we reached a new point
            if (valueMode) {
                obj.setOffset(obj.getOffset() + data::QpNum(valueStr) * sign);
                // new default sign '+'
                sign = 1;
                // new default value 1
                valueStr = "1";
                valueMode = false;
                varNameMode = false;
            } else if (varNameMode) {
                
                if (scientificMode){
                    AuxSciExp=exponentStr;
                    CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_COEF);
                    if(AuxSciCoef==0 && AuxSciExp==0) valueStr="0";
                    AuxSciCoef=0;
                    AuxSciExp=0;
                    if(exponentStr.size()==1 && (exponentStr[0]=='+' || exponentStr[0]=='-'))
                        throw utils::ParserException("Unsupported scientific notation. Did you use e/E as variable name? " + constraint);
                    exponentStr="";
                    scientificMode = false;
                }
                else CheckDigits(valueStr);
                
                val = valueStr;
                obj.setObjElement(getVarByName(varName, vars).getIndex(), val * sign);
                // new default sign '+'
                sign = 1;
                // new default value 1
                valueStr = "1";
                valueMode = false;
                varNameMode = false;
            }
            // the new sign
            c == '+' ? sign = 1 : sign = -1;
            valueStr = "1";
        } else if (isdigit(c) || c == '.') {
            if (valueMode) {
                valueStr += c;
                if (scientificMode){
                    if(c == '.')
                        throw utils::ParserException("only use integers as exponents for scientific notation " + constraint);
                    exponentStr += c;
                }
            } else if (varNameMode) {
                varName += c;
            } else {
                valueMode = true;
                valueStr = c;
            }
        } else /*if ((c) || c == '_')*/{
            if (varNameMode) {
                varName += c;
            } else if (valueMode || scientificMode) {
		if((c == 'e' || c == 'E')
                    && i+1 < objectiveFunction.size() && (objectiveFunction[i+1] == '-' || objectiveFunction[i+1] =='+' || isdigit(objectiveFunction[i+1]))){
                    bool PossibleVariable=false;
                    if(isdigit(objectiveFunction[i+1])){
                        //Maybe this is in fact a variable name like E1 or e2411; Note that variable names like e1d are forbidden
                        PossibleVariable=true;
                        int loc=i+1;
                        while(loc < objectiveFunction.size() && objectiveFunction[loc]!= '+' && objectiveFunction[loc]!= '-' && objectiveFunction[loc]!= ' ' && objectiveFunction[loc]!= '<' && objectiveFunction[loc]!= '>' && objectiveFunction[loc]!= '='){
                            if(isdigit(objectiveFunction[loc])) loc++;
                            else{
                                PossibleVariable=false;
                                break;
                            }
                        }
                    }
                    if(!PossibleVariable){
                        if((c == 'e' || c == 'E') && ( i+1 < objectiveFunction.size() && (objectiveFunction[i+1] == '-' || objectiveFunction[i+1] =='+' || isdigit(objectiveFunction[i+1])))){
                            if (scientificMode){
                                throw utils::ParserException("Unsupported Double Exponent in scientific notation. Also do not use e/E as variable name!");
                            }
                            CheckDigits(valueStr,1);
                            AuxSciCoef=valueStr;
                            scientificMode=true;
                            valueStr+=c;
                            if(points[i+1]=='-' || points[i+1]=='+'){
                                exponentStr=points[i+1];
                                valueStr+=points[i+1];
                                i++;    //skip explicitly reading the sign of exponent; but add it to the strings to be handeled by QpRHS/QpVal
                                continue;
                            }
                        }
                    }
                    else{//Found possible variable like e121 or E7
                        varName = c;
                        valueMode = false;
                        varNameMode = true;
                    }
            	}
            	else{
                    // end value mode and start verNameMode
                    varName = c;
                    valueMode = false;
                    scientificMode = false;
                    varNameMode = true;
                }
            } else {
                varName = c;
                varNameMode = true;
            }
        }
        if (i == objectiveFunction.size() - 1) {
            if (scientificMode){
                AuxSciExp=exponentStr;
                CheckScientificNotation(AuxSciCoef,AuxSciExp,ST_COEF);
                if(AuxSciCoef==0 && AuxSciExp==0) valueStr="0";
                AuxSciCoef=0;
                AuxSciExp=0;
                if(exponentStr.size()==1 && (exponentStr[0]=='+' || exponentStr[0]=='-'))
                    throw utils::ParserException("Unsupported scientific notation. Did you use e/E as variable name? " + constraint);
                exponentStr="";
                scientificMode = false;
            }
            else CheckDigits(valueStr);
            val = valueStr;
            
            
            if(valueMode) obj.setOffset(obj.getOffset() + data::QpNum(valueStr) * sign);
            else obj.setObjElement(this->getVarByName(varName, vars).getIndex(), val * sign);
            sign = 1;
            valueStr = "1";
            valueMode = false;
            varNameMode = false;
        }
    }
}
}
