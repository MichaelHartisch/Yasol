/*
*
* Solver: Log.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/Log.hpp"
#include "Settings/Settings.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

namespace utils {

std::list<Logger*> Logger::loggerList;

Logger::~Logger() {
	log(LOG_DEBUG, "Logger::~Logger", "Logger " + name + ", signing off");
	if (allocatedStream)
		delete target;
}

std::string stringLevel(level l) {
	switch (l) {
	case LOG_ERROR:
		return "ERROR ";
	case LOG_INFO:
		return "INFO  ";
	case LOG_DEBUG:
		return "DEBUG ";
	case LOG_INSANE:
		return "INSANE";
	default:
		return "UNKNOWN LOG LEVEL";
	}
}

level levelString(std::string level) {
	if (level == std::string("ERROR"))
		return LOG_ERROR;
	else if (level == std::string("INFO"))
		return LOG_INFO;
	else if (level == std::string("DEBUG"))
		return LOG_DEBUG;
	else if (level == std::string("INSANE"))
		return LOG_INSANE;
	else
		Logger::globalLog(LOG_ERROR, "levelString",
				"unknown log level. returning INSANE");
	return LOG_INSANE;
}

void Logger::newLogger(std::ostream* target, level debugLevel, std::string name) {
	Logger* logger = new Logger(target, debugLevel, name);
	loggerList.push_back(logger);
}

void Logger::newLogger(std::string filename, level debugLevel, std::string name) {
	Logger* logger = new Logger(filename, debugLevel, name);
	loggerList.push_back(logger);
}

void Logger::addTimeCode(std::string loggerName) {
	std::list<Logger*>::iterator it;
	for (it = loggerList.begin(); it != loggerList.end(); ++it)
		if ((*it)->getName() == loggerName)
			(*it)->timeCode = true;
}

void Logger::remove(std::string name) {
	std::list<Logger*>::iterator it;
	for (it = loggerList.begin(); it != loggerList.end(); ++it)
		if ((*it)->getName() == name) {
			Logger* eintagsfliege = *it;
			loggerList.erase(it);
			delete eintagsfliege;
			it = loggerList.begin();
		}
}

void Logger::cleanup() {
	globalLog(LOG_DEBUG, "Logger::cleanup", "Logger Apocalypse Now!");
	std::list<Logger*>::iterator it;
	for (it = loggerList.begin(); it != loggerList.end(); it
			= loggerList.begin()) {
		Logger* tmpLog = *it;
		loggerList.erase(it);
		delete tmpLog;
	}
}

void Logger::globalLog(level debugLevel, std::string source, std::string str) {
	std::list<Logger*>::iterator it;
	for (it = loggerList.begin(); it != loggerList.end(); it++)
		(*it)->log(debugLevel, source, str);
}

void Logger::log(level debugLevel, std::string str) {
	log(debugLevel, "unknown origin", str);
}

void Logger::log(level debugLevel, std::string source, std::string str) {
	if (debugLevel <= logLevel) {
		if (timeCode) {
			time_t now = time(0);
			char timestamp[22];
			strftime(timestamp, 22, "%d.%m.%Y - %H:%M:%S", localtime(&now));
			(*target) << timestamp << " ";
		}
		if(TABLE_MODE){
			(*target) << str << std::endl;
		}else{
			(*target) << "<" << name << "> (" << stringLevel(debugLevel) << ") ["
				<< source << "] " << str << std::endl;
		}
	}
}

std::string Logger::getName() {
	return name;
}

Logger::Logger(std::ostream* target_, level debugLevel_, std::string name_) :
	target(target_), logLevel(debugLevel_), timeCode(false), name(name_),
			allocatedStream(false) {
}

Logger::Logger(std::string filename_, level debugLevel_, std::string name_) :
	logLevel(debugLevel_), timeCode(false), name(name_), allocatedStream(true) {
	target = new std::ofstream(filename_.c_str(), std::ios::out
			| std::ios::trunc);
}

}
