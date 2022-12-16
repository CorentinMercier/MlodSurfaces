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

#ifndef PN_H
#define PN_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <fstream>
#include <cstring>
#include "timer.h"
#include <cmath>

#include <immintrin.h>

#include "../dual-contouring/dual_contouring.h"
#include "../../lib/tinyply.h"
#include "Fastapss.h"

//#define CHECK

using namespace std;

///
/// \brief The Pn class
/// Class handling the pointsets
///

class Pn
{
public:
	Pn();
	explicit Pn(std::vector<glm::vec3> points) { initialize(points); }
	explicit Pn(std::vector<glm::vec3> points, std::vector<glm::vec3> normals, glm::vec3 offset = glm::vec3(0), float norme = 1.f) { initialize(points); m_normals = std::move(normals); m_offset = offset; m_norme = norme;}
	virtual ~Pn(){m_positions.clear(); m_normals.clear(); m_density.clear();}
	Pn(const Pn &p);
	Pn(Pn&& p) = default;
	Pn &operator =(Pn &p);
	void clear();
	void load(const string &filename, float factor, bool do_normalize = true, glm::vec3 offset = glm::vec3(0), float norme = 1.f);
	void save(const string & filename) const;
	double initialize(unsigned int nbOfElements, glm::vec3 mini, glm::vec3 maxi, unsigned sizeOfGrid, bool orderWithMorton = true);
	double initialize(unsigned int nbOfElements, Pn &pointSetToFit, float factor, unsigned sizeOfGrid, bool orderWithMorton = true);
	void initialize(std::vector<glm::vec3> points);
	void resize(unsigned int nbOfElements) {
		m_positions.resize(nbOfElements);
		m_normals.resize(nbOfElements);
	}
	void normalize();
	void computeMinMax();
	void computeMeanGap(float factor);

	void add(const glm::vec3 pt, const glm::vec3 n);
	void add(const vector<glm::vec3> pt, const vector<glm::vec3> n);

	const glm::vec3 &operator[](const unsigned int element) const;
	glm::vec3 &operator[](const unsigned int element);
	glm::vec3 &normal(const unsigned int element);
	float density(const unsigned int element);

	glm::vec3 *positionData();
	glm::vec3 *normalData();

	const vector<glm::vec3>& getPositions() const {return m_positions;}
	const vector<glm::vec3>& getNormals() const {return m_normals;}

	double order(unsigned sizeOfGrid = 64);
	void computeMortonOrder(unsigned sizeOfGrid , std::vector<std::pair<uint32_t, unsigned>> & mortonCode);
	void orderWithColor(unsigned nbPts, Pn &newSkeleton, vector<glm::vec3> & colors, float factor);
	void clustering(vector<glm::vec3> & colors);

	unsigned int size() const;
	unsigned int validSize() const {return m_validSize;}
	void setValidSize(unsigned int newValidSize){m_validSize = newValidSize;}

	void resample(unsigned int nbOfPoints);

	void upsample(unsigned int divisionFactor, Pn &skeleton);

	glm::vec3 getMini() const {return m_mini;}
	glm::vec3 getMaxi() const {return m_maxi;}
	glm::vec3 getOffset()const {return m_offset;}
	float getNorme() const {return m_norme;}

	void computePointsDensity();

	static glm::vec3 GetColour(float v,float vmin,float vmax);

private:

	vector<glm::vec3> m_positions = vector<glm::vec3>(0);
	vector<glm::vec3> m_normals = vector<glm::vec3>(0);

	vector<float> m_density = vector<float>(0);

	glm::vec3 m_mini;
	glm::vec3 m_maxi;

	float m_norme = 1.f;
	glm::vec3 m_offset = glm::vec3(0, 0, 0);

	float m_meanGap;
	unsigned int m_validSize;
};

#endif // PN_H
