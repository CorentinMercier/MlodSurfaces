// --------------------------------------------------------------------------
// Source code provided FOR REVIEW ONLY, as part of the submission entitled
// "Moving Level-of-Detail Surfaces".
//
// A proper version of this code will be released if the paper is accepted
// with the proper licence, documentation and bug fix.
// Currently, this material has to be considered confidential and shall not
// be used outside the review process.
//
// All right reserved. The Authors
// --------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <chrono>

///
/// \brief The Timer class Class to improve timer readability
///

class CPUTimer
{
	public:
	CPUTimer() { start(); }

	inline void start() { _start = std::chrono::high_resolution_clock::now(); }
	inline void reset() { _start = std::chrono::high_resolution_clock::now(); }
	inline void stop() { _stop = std::chrono::high_resolution_clock::now(); }
	inline double elapsed() { stop();  return std::chrono::duration<double>(_stop - _start).count(); }

	inline double printElapsed()
	{
		double delta = elapsed();
		std::cout << delta << " sec" << std::endl;
		return delta;
	}

	inline double printElapsed(const std::string & s)
	{
		double delta = elapsed() * 1000;
		std::cout << s.c_str() << ": " << delta << " ms" << std::endl;
		return delta;
	}

	inline double printElapsedAndReset()
	{
		double delta = printElapsed();
		reset();
		return delta;
	}

	inline double printElapsedAndReset(const std::string & s)
	{
		double delta = printElapsed(s);
		reset();
		return delta;
	}

	private:
	std::chrono::time_point<std::chrono::high_resolution_clock> _start, _stop;
};
