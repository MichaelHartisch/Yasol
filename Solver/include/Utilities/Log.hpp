/*
*
* Solver: Log.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef LOG_H_
#define LOG_H_

#include <string>
#include <list>

namespace utils {

typedef enum LEVEL
{
	LOG_ERROR, LOG_INFO, LOG_DEBUG, LOG_INSANE
} level;

/** Printable Log Levels.
 * This function returns a human-readable string representing the
 * given log level.
 * @param l enum type
 * @return string representation of enum type
 */
std::string stringLevel(LEVEL l);

/** Log Level String Interpreter.
 * This function returns the log level corresponding
 * to the given string representation. It is used to
 * read log levels from a human-readable ConfigFile.
 * @param level log level string
 * @return enumeration type
 */
level levelString(std::string level);

/** Log Class.
 * This class implements a logging system. It provides a registry
 * for a set of loggers, functions for adding and deleting
 * loggers and for logging to all currently registered loggers.
 */
class Logger {
public:
	~Logger();

	/** Create a new Stream Logger.
	 * This method creates a new logger for the given stream and registers it.
	 * Only log messages with a log level lower than or equal to debugLevel
	 * will be written to the logger. The logger will be identified by name.
	 * @param target outstream of logger
	 * @param debugLevel maximal log level that gets logged by this logger
	 * @param name name of the logger
	 */
	static void newLogger(std::ostream* target, level debugLevel,
			std::string name);

	/** Create a new File Logger.
	 * This method creates a new logger for the given file and registers it.
	 * Only log messages with a log level lower than or equal to debugLevel
	 * will be written to the logger. The logger will be identified by name.
	 * @param filename path of filename logger writes to
	 * @param debugLevel maximal log level that gets logged by this logger
	 * @param name name of the logger
	 */
	static void newLogger(std::string filename, level debugLevel,
			std::string name);

	/** Add a time code to a logger.
	 * This method sets the timecode flag for the logger with the
	 * given name.
	 * @param loggerName name of the logger
	 */
	static void addTimeCode(std::string loggerName);

	/** Add a thread field to a logger.
	 * This method sets the threadField flag for the logger with the
	 * given name.
	 * @param loggerName name of the logger
	 */
	static void addThreadField(std::string loggerName);

	/**
	 * Prints hostname of the machine the current process is running to the logger.
	 * @param loggerName name of the logger
	 */
	static void addIP(std::string loggerName);

	/**
	 * Sets the hostname of the machine the current process is running.
	 * @param hostname string of hostname
	 */
	void setHostname(std::string hostname);

	/** Remove Logger.
	 * This method unregisters and destroys the logger with the given name.
	 * @param name name of the logger
	 */
	static void remove(std::string name);

	/** Remove all loggers.
	 * This method unregisters and destroys all currently registered loggers.
	 */
	static void cleanup();

	/** Log globally.
	 * This method sends the given log message to all registered loggers.
	 * The namespace and method in which this method is called should be
	 * written into the source argument.
	 * @param debugLevel log level
	 * @param source method where log call originates
	 * @param str message to log
	 */
	static void
			globalLog(level debugLevel, std::string source, std::string str);

	/** Get the Name of this Logger.
	 */
	std::string getName();

	/*
	static std::string NumberToString(double number) {

		std::stringstream NumberString;
		NumberString << number;
		std::string Number = NumberString.str();

		return Number;
	}*/

private:
	/** Private Constructor.
	 * This constructor corresponds to the newLogger method. It is private to prevent
	 * the creation of loggers outside this class.
	 * @param target_ outstream
	 * @param debugLevel_ maximal log level
	 * @param name_ name of the logger
	 */
	Logger(std::ostream* target_, level debugLevel_, std::string name_);

	/** Private Constructor.
	 * This constructor corresponds to the newLogger method. It is private to prevent
	 * the creation of loggers outside this class.
	 * @param filename_ path to logger file name
	 * @param debugLevel_ maximal log level
	 * @param name_ name of the logger
	 */
	Logger(std::string filename_, level debugLevel_, std::string name_);

	/* Log locally.
	 * This method sends the given log message to the stream connected
	 * to this logger.
	 * If possible, this message should not be used because it does not
	 * print out information about where the log message is coming from.
	 * @param debugLevel maximal log level
	 * @param str message to log
	 */
	void log(level debugLevel, std::string str);

	/* Log locally.
	 * This method sends the given log message to the stream connected
	 * to this logger.
	 * The namespace and method in which this method is called should be
	 * written into be source argument.
	 * @param debugLevel log level
	 * @param source method where this log call originates
	 * @param str message to log
	 */
	void log(level debugLevel, std::string source, std::string str);

	/** Check for name collision.
	 * This checks whether a logger of the given name already exists.
	 * @param name name of logger
	 * @return exists the given logger?
	 */
	static bool nameExistence(std::string name);

	/** Target stream.
	 * This is the stream this logger uses for output.
	 */
	std::ostream* target;

	/** Log level.
	 * All messages sent to this logger with a log level higher than
	 * logLevel will be ignored.
	 */
	level logLevel;

	/** Time Code Flag.
	 * If set, the logger will add a time code in front of all log messages.
	 */
	bool timeCode;

	/** The name of this logger.
	 */
	std::string name;

	/** Allocation Indicator.
	 * Indicates whether the ostream was self-allocated. If so, the memory
	 * needs to be manually freed upon destruction of this logger.
	 */
	bool allocatedStream;

	/** Logger registry.
	 * This list contains all loggers that are currently active.
	 */
	static std::list<Logger*> loggerList;

};
}

#endif /* LOG_H_ */
