// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    Moving Level-of-Detail Surfaces.
//    C. Mercier, T. Lescoat, P. Roussillon, T. Boubekeur, and J-M. Thiery
//    ACM Transaction On Graphics 2022
//    DOI: 10.1145/3528223.3530151
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.
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
