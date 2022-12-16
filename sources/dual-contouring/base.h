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

#include <Eigen/Dense>

#ifdef WIN32
#include <intrin.h>
#undef near
#undef far
inline void debug_break() { __debugbreak(); }
#else
inline void debug_break() {}
#endif

#define ASSERT( condition, message ) \
	if(!(condition)) { printf("ASSERT: %s\n", message); debug_break(); }


template<typename Scalar>
Eigen::Matrix<Scalar, 4, 4> translation4x4(const Eigen::Matrix<Scalar, 3, 1>& v)
{
	Eigen::Matrix<Scalar, 4, 4> mat = Eigen::Matrix<Scalar, 4, 4>::Identity();
	for(int i = 0; i < 3; i++)
		mat(i, 3) = v(i);
	return mat;
}

template<typename Scalar>
Eigen::Matrix<Scalar, 4, 4> scale4x4(Scalar v)
{
	Eigen::Matrix<Scalar, 4, 4> mat = Eigen::Matrix<Scalar, 4, 4>::Identity() * v;
	mat(3, 3) = Scalar(1);
	return mat;
}


template<typename Scalar>
Eigen::Matrix<Scalar, 4, 1> to_hpos(const Eigen::Matrix<Scalar, 3, 1>& v)
{
	return { v[0], v[1], v[2], Scalar(1) };
}

template<typename Scalar>
Eigen::Matrix<Scalar, 3, 1> from_hpos(const Eigen::Matrix<Scalar, 4, 1>& v)
{
	return { v[0] / v[3], v[1] / v[3], v[2] / v[3] };
}

inline Eigen::Vector3f color_ramp(float v)
{
	if(v > 1) v = 1;
	if(v < -1) v = -1;
	auto sat = [](float t) -> float { return t > 1.0f ? 1.0f : t < 0.0f ? 0.0f : t; };
	return { sat(1.5f - fabs(2 * v - 1)), sat(1.5f - fabs(2 * v)), sat(1.5f - fabs(2 * v + 1)) };
}


template<typename I>
struct integer_range
{
	struct integer_iterator
	{
		I t;

		I operator*() const { return t; }
		integer_iterator& operator++() { t++; return *this; }
		integer_iterator operator++(int) { return ++integer_iterator{ t }; }
		bool operator==(const integer_iterator& other) const { return t == other.t; }
		bool operator!=(const integer_iterator& other) const { return t != other.t; }
	};

	I begin_value, end_value;
	integer_range(I begin, I end) : begin_value(std::min(begin, end + 1)), end_value(begin != end ? std::max(end, begin + 1) : end) {}
	integer_iterator begin() const { return { begin_value }; }
	integer_iterator end() const { return { end_value }; }
};

template<typename T>
struct pointer_range
{
	T* begin_ptr;
	T* end_ptr;

	pointer_range(T* begin, T* end) : begin_ptr(begin), end_ptr(end) {}
	T* begin() const { return begin_ptr; }
	T* end() const { return end_ptr; }
};

inline bool vec3i_less_xyz(const Eigen::Vector3i& a, const Eigen::Vector3i& b)
{
	for(int i = 0; i < 3; i++)
		if(a[i] != b[i])
			return a[i] < b[i];
	return false;
}

template<typename Range>
inline bool is_empty(const Range& r) { return r.begin() == r.end(); }

enum class direction { X = 0, Y = 1, Z = 2 };
