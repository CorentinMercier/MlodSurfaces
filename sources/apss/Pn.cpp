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

#include "Pn.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <set>
#include <random>


glm::vec3 Pn::GetColour(float v,float vmin,float vmax)
{
   glm::vec3 c = glm::vec3(1.f,1.f,1.f); // white
   float dv;
   if (v < vmin)
	  v = vmin;
   if (v > vmax)
	  v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
	  c.x = 0;
	  c.y = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
	  c.x = 0;
	  c.z = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
	  c.x = 4 * (v - vmin - 0.5 * dv) / dv;
	  c.z = 0;
   } else {
	  c.y = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
	  c.z = 0;
   }

   return(c);
}

Pn::Pn()
{

}

Pn::Pn(const Pn &p)
{
	m_positions.resize(p.size());
	m_normals.resize(p.size());
#pragma omp parallel for
	for (int i=0; i<(int)m_positions.size(); i++)
		m_positions[i] = p.m_positions[i];
#pragma omp parallel for
	for (int i=0; i<(int)p.getNormals().size(); i++)
		m_normals[i] = p.m_normals[i];
	m_mini = p.m_mini;
	m_maxi = p.m_maxi;
	m_meanGap = p.m_meanGap;
	m_offset = p.m_offset;
	m_norme = p.m_norme;
	m_validSize = p.m_validSize;
}

Pn& Pn::operator =(Pn &p)
{
	m_positions.resize(p.size());
	m_normals.resize(p.size());
#pragma omp parallel for
	for (int i=0; i<(int)m_positions.size(); i++)
		m_positions[i] = p[i];
#pragma omp parallel for
	for (int i=0; i<(int)p.getNormals().size(); i++)
		m_normals[i] = p.normal(i);
	m_mini = p.m_mini;
	m_maxi = p.m_maxi;
	m_meanGap = p.m_meanGap;
	m_offset = p.m_offset;
	m_norme = p.m_norme;
	m_validSize = p.m_validSize;
	return *this;
}

void Pn::clear()
{
	m_positions.clear();
	m_validSize = 0;
	m_normals.clear();
	m_density.clear();
}

bool hasEnding(const string & path, const string & extension)
{
	if (path.length()<extension.length()) return false;
	return (path.compare(path.length()-extension.length(), extension.length(), extension) == 0);
}

void Pn::load(const string & filename, float factor, bool do_normalize, glm::vec3 offset, float norme)
{
	if (hasEnding(filename, "pn"))
	{
		unsigned int surfelSize = 6;
		FILE * in = fopen (filename.c_str (), "rb");
		if (in == nullptr) {
			cout << filename << " is not a valid PN file." << endl;
			return;
		}
		size_t READ_BUFFER_SIZE = 1000; // for example...
		float * pn = new float[surfelSize*READ_BUFFER_SIZE];
		m_positions.clear ();
		m_normals.clear ();
		while (!feof (in)) {
			unsigned numOfPoints = fread (pn, 4, surfelSize*READ_BUFFER_SIZE, in);
			for (unsigned int i = 0; i < numOfPoints; i += surfelSize) {
				m_positions.push_back (glm::vec3 (pn[i], pn[i+1], pn[i+2]));
				m_normals.push_back (glm::vec3 (pn[i+3], pn[i+4], pn[i+5]));
			}

			if (numOfPoints < surfelSize*READ_BUFFER_SIZE) break;
		}
		delete [] pn;
		fclose (in);
	}
	else if (hasEnding(filename, "ply"))
	{
		std::ifstream ss(filename, std::ios::binary);
		if (ss.fail()) {
			cout << filename << " is not a valid ply file." << endl;
			return;
		}
		tinyply::PlyFile pf;
		pf.parse_header(ss);
		std::shared_ptr<tinyply::PlyData> vertices, normals;
		try { vertices = pf.request_properties_from_element("vertex", { "x", "y", "z" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
		bool thereAreNormals = true;
		try { normals = pf.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; thereAreNormals = false;}
		pf.read(ss);
		m_positions.resize(vertices->count);
		std::memcpy(m_positions.data(), vertices->buffer.get(), vertices->buffer.size_bytes());
		if (thereAreNormals)
		{
			m_normals.resize(normals->count);
			std::memcpy(m_normals.data(), normals->buffer.get(), normals->buffer.size_bytes());
		}
		ss.close();
	}
	if (do_normalize)
		normalize();
	else
	{
		for (unsigned i=0; i<size(); i++)
			m_positions[i] = (m_positions[i] - offset) / norme;
		m_offset = offset;
		m_norme = norme;
	}
	computeMeanGap(factor);
	m_validSize = m_positions.size();
}

void Pn::computeMeanGap(float factor)
{
	cout << "Factor: " << factor << endl;
	float r =0.5f;
	float surface = 4 * M_PI * r * r;
	m_meanGap = (cbrt(surface / (float(m_positions.size()))) * 2.f) * factor;
}

void Pn::add(const glm::vec3 pt, const glm::vec3 n)
{
	for (unsigned int i=0; i<m_positions.size(); i++)
	{
		if (pt == m_positions[i] && n == m_normals[i])
			return;
	}
	m_positions.push_back(pt);
	m_validSize = m_positions.size();
	m_normals.push_back(n);
}

void Pn::add(const vector<glm::vec3> pt, const vector<glm::vec3> n)
{
	if (pt.size() != n.size())
	{
		cerr << "Error, different number of points and normals" << endl;
		return;
	}
	for (unsigned i=0; i<pt.size(); i++)
	{
		m_positions.push_back(pt[i]);
		m_normals.push_back(n[i]);
	}
	m_validSize = m_positions.size();
}

void Pn::save(const string & filename) const
{
	if (hasEnding(filename, "pn"))
	{
		FILE * outfile = fopen (filename.c_str (), "wb");
		if (outfile == NULL) {
			cout << filename << " is not a valid PN file." << endl;
			return;
		}
		for(unsigned int pIt = 0 ; pIt < m_positions.size() ; ++pIt) {
			glm::vec3 realPosition = (m_positions[pIt] * m_norme) + m_offset;
			fwrite (&(realPosition) , sizeof(float), 3, outfile);
			fwrite (&(m_normals[pIt]) , sizeof(float), 3, outfile);
	}
	fclose (outfile);
	}
	else if (hasEnding(filename, "ply"))
	{
		fstream file;
		file.open(filename, ios_base::out | ios_base::binary);
		if (!file.is_open())
		{
			std::cerr << "Impossible to write file " << filename << std::endl;
			return;
		}
		file << "ply" << endl;
		file << "format binary_little_endian 1.0" << endl;
		file << "element vertex " << int(m_positions.size()) << endl;
		file << "property float x" << endl;
		file << "property float y" << endl;
		file << "property float z" << endl;
		if (m_normals.size() != 0)
		{
			file << "property float nx" << endl;
			file << "property float ny" << endl;
			file << "property float nz" << endl;
		}
		file << "end_header" << endl;
		for (unsigned int i=0; i<m_positions.size(); i++)
		{
			glm::vec3 realPosition = (m_positions[i] * m_norme) + m_offset;
			file.write(reinterpret_cast<const char *>(&(realPosition)), sizeof (glm::vec3));
			if (m_normals.size() != 0)
				file.write(reinterpret_cast<const char *>(&(m_normals[i])), sizeof (glm::vec3));
		}
		file.close();
	}
}

void Pn::initialize(std::vector<glm::vec3> points)
{
	m_positions = std::move(points);
	m_validSize = m_positions.size();
	m_normals.clear();
	m_density.clear();
	m_norme = 1;
	m_offset = { 0, 0, 0 };
	computeMinMax();
	computeMeanGap(1.0f);
}

void Pn::normalize()
{
	m_offset = m_positions[0];
	m_maxi = m_positions[0];
	m_mini = m_positions[0];
	for (unsigned i=1; i<m_positions.size(); i++)
	{
		m_offset += m_positions[i];
		m_maxi.x = max(m_positions[i].x, m_maxi.x);
		m_maxi.y = max(m_positions[i].y, m_maxi.y);
		m_maxi.z = max(m_positions[i].z, m_maxi.z);
		m_mini.x = min(m_positions[i].x, m_mini.x);
		m_mini.y = min(m_positions[i].y, m_mini.y);
		m_mini.z = min(m_positions[i].z, m_mini.z);
	}
	m_offset /= m_positions.size();
	m_maxi -= m_offset;
	m_mini -= m_offset;

	const float smallFactor = glm::length(m_maxi - m_mini) * 0.005f;
	m_maxi += glm::vec3(smallFactor);
	m_mini -= glm::vec3(smallFactor);

	m_norme = max(max(max(max(max(fabs(m_maxi.x), fabs(m_maxi.y)), fabs(m_maxi.z)), fabs(m_mini.x)), fabs(m_mini.y)), fabs(m_mini.z));

#pragma omp parallel for
	for (int i=0; i<(int)m_positions.size(); i++)
		m_positions[i] = (m_positions[i] - m_offset) / m_norme;
}

void Pn::computeMinMax()
{
	m_mini = m_positions[0];
	m_maxi = m_positions[0];
	for (unsigned int i=1; i<m_positions.size(); i++)
	{
		for (unsigned int j=0; j<3; j++)
		{
			if (m_positions[i][j] < m_mini[j])
				m_mini[j] = m_positions[i][j];
			if (m_positions[i][j] > m_maxi[j])
				m_maxi[j] = m_positions[i][j];
		}
	}
}

double Pn::initialize(unsigned int nbOfElements, glm::vec3 mini, glm::vec3 maxi, unsigned sizeOfGrid, bool orderWithMorton)
{
	m_positions.resize(nbOfElements);
	m_validSize = nbOfElements;
#pragma omp parallel for
	for (int i=0; i<(int)nbOfElements; i++)
		m_positions[i] = glm::vec3(
					mini.x + (maxi.x - mini.x) * (double)(rand())/(double)(RAND_MAX),
					mini.y + (maxi.y - mini.y) * (double)(rand())/(double)(RAND_MAX),
					mini.z + (maxi.z - mini.z) * (double)(rand())/(double)(RAND_MAX)
					);
	m_mini = mini;
	m_maxi = maxi;
	return (orderWithMorton) ? order(sizeOfGrid) : 0.0;
}

double Pn::initialize(unsigned int nbOfElements, Pn & pointSetToFit, float factor, unsigned sizeOfGrid, bool orderWithMorton)
{
	m_meanGap = pointSetToFit.m_meanGap;
	m_offset = pointSetToFit.m_offset;
	m_norme = pointSetToFit.m_norme;

	//Better fitted
	m_positions.resize(nbOfElements);
	m_validSize = nbOfElements;
	m_normals.resize(nbOfElements);
	computeMeanGap(factor);
	unsigned int size = pointSetToFit.size();
	cout << "Size of pointset to fit: " << size << endl;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<unsigned int> dist(0, size - 1);
	for (unsigned int i=0; i<nbOfElements; i++)
	{
		unsigned int randomPoint = dist(mt);
		m_positions[i] = pointSetToFit[randomPoint];
		m_normals[i] = pointSetToFit.normal(randomPoint);
		glm::vec3 e1;
		if (glm::dot(pointSetToFit.normal(randomPoint), glm::vec3(0, 0, 1)) != 0)
			e1 = glm::normalize(glm::cross(pointSetToFit.normal(randomPoint), glm::vec3(0, 0, 1)));
		else
			e1 = glm::normalize(glm::cross(pointSetToFit.normal(randomPoint), glm::vec3(0, 1, 0)));
		glm::vec3 e2 = glm::normalize(glm::cross(pointSetToFit.normal(randomPoint), e1));
		float randR = m_meanGap * sqrt(rand() / (float)RAND_MAX);
		float randTheta = rand() / float(RAND_MAX) * 2 * M_PI;
		glm::vec3 randomVector = randR * cos(randTheta) * e1 + randR * sin(randTheta) * e2;
		m_positions[i] += randomVector;
	}
	computeMinMax();
	return (orderWithMorton) ? order(sizeOfGrid) : 0.0;
}


void Pn::computeMortonOrder(unsigned sizeOfGrid , std::vector<std::pair<uint32_t, unsigned>> & mortonCode) {
	mortonCode.resize(size());

	computeMinMax();

	float mini = min(min(m_mini.x, m_mini.y), m_mini.z);
	float maxi = max(max(m_maxi.x, m_maxi.y), m_maxi.z);
	float norme = maxi - mini;

	//#pragma omp parallel for
		for (unsigned int i=0; i<size(); i++)
		{
			mortonCode[i].first = 0;
			mortonCode[i].second = i;
			uint32_t mortonValue1 = 0, mortonValue2 = 0, mortonValue3 = 0;
			mortonValue1 = (m_positions[i].x - mini) / norme * sizeOfGrid;
			mortonValue2 = (m_positions[i].y - mini) / norme * sizeOfGrid;
			mortonValue3 = (m_positions[i].z - mini) / norme * sizeOfGrid;

			mortonCode[i].first = _pdep_u32(mortonValue1, 0x92492492) | _pdep_u32(mortonValue2, 0x49249249) | _pdep_u32(mortonValue3, 0x24924924);

		}
		std::sort(mortonCode.begin(), mortonCode.end(), [](std::pair<uint32_t, unsigned> val1, std::pair<uint32_t, unsigned> val2) {
			return val1.first < val2.first;
		});
}

double Pn::order(unsigned sizeOfGrid)
{
	CPUTimer t;
	t.start();

	std::vector<unsigned> nbOfPtsPerCluster(sizeOfGrid * sizeOfGrid * sizeOfGrid, 0);

	//With Morton code
	std::vector<std::pair<uint32_t, unsigned>> mortonCode(size());

	std::set< uint32_t > codes;

	float mini = min(min(m_mini.x, m_mini.y), m_mini.z);
	float maxi = max(max(m_maxi.x, m_maxi.y), m_maxi.z);
	float norme = maxi - mini;

	for (unsigned int i=0; i<size(); i++)
	{
		mortonCode[i].first = 0;
		mortonCode[i].second = i;
		const uint32_t mortonValue1 = glm::clamp<uint32_t>(std::floor((m_positions[i].x - mini) / norme * sizeOfGrid), 0, sizeOfGrid - 1);
		const uint32_t mortonValue2 = glm::clamp<uint32_t>(std::floor((m_positions[i].y - mini) / norme * sizeOfGrid), 0, sizeOfGrid - 1);
		const uint32_t mortonValue3 = glm::clamp<uint32_t>(std::floor((m_positions[i].z - mini) / norme * sizeOfGrid), 0, sizeOfGrid - 1);

		nbOfPtsPerCluster[mortonValue1 + mortonValue2 * sizeOfGrid + mortonValue3 * sizeOfGrid * sizeOfGrid]++;

		mortonCode[i].first = _pdep_u32(mortonValue1, 0x92492492) | _pdep_u32(mortonValue2, 0x49249249) | _pdep_u32(mortonValue3, 0x24924924);

		codes.insert( mortonCode[i].first );
	}
	std::sort(mortonCode.begin(), mortonCode.end(), [](std::pair<uint32_t, unsigned> val1, std::pair<uint32_t, unsigned> val2) {
		return val1.first < val2.first;
	});
	std::vector<glm::vec3> newPositions(size());
	for (unsigned int i=0; i<size(); i++)
		newPositions[i] = m_positions[mortonCode[i].second];
	for (unsigned int i=0; i<size(); i++)
		m_positions[i] = newPositions[i];

	double elapsed = t.printElapsed("Morton code computed in : ");
	std::cout << " ...  (" << codes.size() << " different codes)" << std::endl;

#ifdef CHECK
	double meanNbOfPts = 0;
	unsigned nbOfCubes = 0;
	for (unsigned int i=0; i<nbOfPtsPerCluster.size(); i++)
	{
		if (nbOfPtsPerCluster[i] != 0)
		{
			meanNbOfPts += nbOfPtsPerCluster[i];
			nbOfCubes++;
		}
	}
	meanNbOfPts /= float(nbOfCubes);
	cout << "grid size : " << sizeOfGrid << " Mean nb of points : " << meanNbOfPts << " grouped into " << nbOfCubes << " morton cubes" << endl;
#endif
	return elapsed;
}

bool intersectEdge(const float dl, const glm::vec3 & n, const glm::vec3 & pt1, const glm::vec3 & pt2, glm::vec3 & intersection)
{
	float denominator = glm::dot(n, pt2 - pt1);
	if (denominator == 0)
		return false;
	float lambda = (dl - glm::dot(n, pt1)) / (denominator);
	if (lambda < 0 || lambda > 1)
	{
		return false;
	}
	intersection = pt1 + lambda * (pt2 - pt1);
	return true;
}

void Pn::orderWithColor(unsigned nbPts, Pn & newSkeleton, vector<glm::vec3> &colors, float factor)
{
	CPUTimer t;
	t.start();
	computeMeanGap(factor);
	//With Morton code
	unsigned power = pow(2, 5);
	std::vector<std::pair<uint32_t, unsigned>> mortonCode(size());
	for (unsigned int i=0; i<size(); i++)
	{
		mortonCode[i].first = 0;
		mortonCode[i].second = i;
		uint32_t mortonValue1 = 0, mortonValue2 = 0, mortonValue3 = 0;
		mortonValue1 = (m_positions[i].x + 1.f) / 2.f * power;
		mortonValue2 = (m_positions[i].y + 1.f) / 2.f * power;
		mortonValue3 = (m_positions[i].z + 1.f) / 2.f * power;

		mortonCode[i].first = _pdep_u32(mortonValue1, 0x92492492) | _pdep_u32(mortonValue2, 0x49249249) | _pdep_u32(mortonValue3, 0x24924924);
	}
	std::sort(mortonCode.begin(), mortonCode.end(), [](std::pair<uint32_t, unsigned> val1, std::pair<uint32_t, unsigned> val2) {
		return val1.first < val2.first;
	});
	unsigned int counter = 0;
	uint32_t previousValue = mortonCode[0].first;
	vector<unsigned> nbPointsPerCluster(1);
	nbPointsPerCluster[0] = 0;
	vector<unsigned int> clusters(size());
	vector<vector<unsigned int> > ptsInClusters;
	vector<unsigned int> ptsInCurrentCluster;
	for (unsigned int i=0; i<size(); i++)
	{
		if (previousValue != mortonCode[i].first)
		{
			counter++;
			nbPointsPerCluster.push_back(0);
			ptsInClusters.push_back(ptsInCurrentCluster);
			ptsInCurrentCluster.clear();
		}
		nbPointsPerCluster[counter]++;
		ptsInCurrentCluster.push_back(mortonCode[i].second);
		previousValue = mortonCode[i].first;
		clusters[i] = counter;
	}
	ptsInClusters.push_back(ptsInCurrentCluster);
	cout << counter+1 << " clusters" << endl;

	t.printElapsedAndReset("Morton : ");

	vector<float> densityPerCluster(counter);
	for (unsigned int i=0; i<counter; i++)
		densityPerCluster[i] = nbPointsPerCluster[i];
	t.printElapsedAndReset("Gaussian : ");

	vector<float> temp(densityPerCluster.begin(), densityPerCluster.end());
	std::nth_element(temp.begin(), temp.begin() + temp.size() / 4, temp.end());
	float firstQuarterDensity = temp[temp.size() / 4];
	std::nth_element(temp.begin(), temp.begin() + temp.size()* 3 / 4, temp.end());
	float thirdQuarterDensity = temp[temp.size() * 3 / 4];

	std::vector<glm::vec3> newPositions(size());
	std::vector<glm::vec3> newNormals(size());
#pragma omp parallel for
	for (int i=0; i<(int)size(); i++)
	{
		newPositions[i] = m_positions[mortonCode[i].second];
		newNormals[i] = m_normals[mortonCode[i].second];
	}
	t.printElapsedAndReset("Assignation of new positions and normals : ");

	//Subsample
	unsigned position = 0;
	vector<unsigned int> possiblePointsToRemove;
	possiblePointsToRemove.reserve(m_positions.size()/4);
	vector<unsigned int> possiblePointsToUpscale;
	possiblePointsToUpscale.reserve(m_positions.size()/4);
	for (unsigned int i=0; i<nbPointsPerCluster.size(); i++)
	{
		if (densityPerCluster[i] > thirdQuarterDensity)
		{
			for (unsigned j=0; j<nbPointsPerCluster[i]; j++)
				possiblePointsToRemove.push_back(position + j);
		}
		else if (densityPerCluster[i] < firstQuarterDensity) {
			for (unsigned j=0; j<nbPointsPerCluster[i]; j++)
				possiblePointsToUpscale.push_back(position + j);
		}
		position += nbPointsPerCluster[i];
	}
	t.printElapsedAndReset("Choosing points : ");
	std::random_shuffle(possiblePointsToRemove.begin(), possiblePointsToRemove.end());
	std::random_shuffle(possiblePointsToUpscale.begin(), possiblePointsToUpscale.end());
	t.printElapsedAndReset("Random shuffle : ");
	unsigned int nbOfPointsToChange = std::min((unsigned)std::min(possiblePointsToRemove.size(), possiblePointsToUpscale.size()), nbPts);
	for (unsigned int i=0; i<nbOfPointsToChange; i++)
	{
		glm::vec3 e1;
		if (glm::dot(newNormals[possiblePointsToUpscale[i]], glm::vec3(0, 0, 1)) != 0)
			e1 = glm::normalize(glm::cross(newNormals[possiblePointsToUpscale[i]], glm::vec3(0, 0, 1)));
		else
			e1 = glm::normalize(glm::cross(newNormals[possiblePointsToUpscale[i]], glm::vec3(0, 1, 0)));
		glm::vec3 e2 = glm::cross(newNormals[possiblePointsToUpscale[i]], e1);
		float randR = m_meanGap * sqrt(rand() / (float)RAND_MAX);
		float randTheta = rand() / float(RAND_MAX) * 2 * M_PI;
		glm::vec3 randomVector = randR * cos(randTheta) * e1 + randR * sin(randTheta) * e2;
		newPositions[possiblePointsToRemove[i]] = newPositions[possiblePointsToUpscale[i]] + randomVector;//(randomVector * sqrt(rand() / (float)RAND_MAX) * m_meanGap); //glm::normalize((newPositions[possiblePointsToUpscale[i]] + randomVector) - newSkeleton[possiblePointsToUpscale[i]]) * r + newSkeleton[possiblePointsToUpscale[i]];
		newNormals[possiblePointsToRemove[i]] = newNormals[possiblePointsToUpscale[i]];
	}
	t.printElapsedAndReset("Change points : ");

#pragma omp parallel for
	for (int i=0; i<(int)size(); i++)
	{
		m_positions[i] = newPositions[i];
		m_normals[i] = newNormals[i];
	}
	t.printElapsed("Resampling time : ");
}

const glm::vec3 &Pn::operator[](const unsigned int element) const
{
	return m_positions[element];
}

glm::vec3 &Pn::operator[](const unsigned int element)
{
	return m_positions[element];
}

glm::vec3 &Pn::normal(const unsigned int element)
{
	return m_normals[element];
}

float Pn::density(const unsigned int element)
{
	return m_density[element];
}

glm::vec3* Pn::positionData()
{
	return m_positions.data();
}

glm::vec3* Pn::normalData()
{
	return m_normals.data();
}

unsigned int Pn::size() const
{
	return (static_cast<unsigned int>(m_positions.size()));
}

void Pn::resample(unsigned int nbOfPoints)
{
	if (nbOfPoints > size())
	{
		std::cerr << "Cannot resample to more than current number of points" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (nbOfPoints == size())
		return;//Already good size
	std::vector<unsigned int> values(size());
	for (unsigned int i=0; i<size(); i++)
		values[i] = i;
	vector<glm::vec3> newPositions, newNormals;
	newPositions.reserve(nbOfPoints);
	newNormals.reserve(nbOfPoints);
	unsigned int currentSize = size();
	for (unsigned int i=0; i<nbOfPoints; i++)
	{
		unsigned int randVal = rand() / (float)RAND_MAX * currentSize;
		newPositions.push_back(m_positions[values[randVal]]);
		newNormals.push_back(m_normals[values[randVal]]);
		values[randVal] = values[currentSize-1];
		currentSize--;
	}
	//Copy new position and normal
	m_positions.resize(nbOfPoints);
	m_validSize = nbOfPoints;
	m_normals.resize(nbOfPoints);
#pragma omp parallel for
	for (int i=0; i<(int)nbOfPoints; i++)
	{
		m_positions[i] = newPositions[i];
		m_normals[i] = newNormals[i];
	}
}

void Pn::upsample(unsigned int divisionFactor, Pn &skeleton)
{
	//Subsample
	unsigned int nbOfPoints = size() / divisionFactor;
	std::vector<unsigned int> values(size());
	for (unsigned int i=0; i<size(); i++)
		values[i] = i;
	vector<glm::vec3> newPositions, newNormals, newSkeleton;
	newPositions.reserve(nbOfPoints);
	newNormals.reserve(nbOfPoints);
	newSkeleton.reserve(nbOfPoints);
	unsigned int currentSize = size();
	for (unsigned int i=0; i<nbOfPoints; i++)
	{
		unsigned int randVal = rand() / (float)RAND_MAX * currentSize;
		newPositions.push_back(m_positions[values[randVal]]);
		newSkeleton.push_back(skeleton[values[randVal]]);
		newNormals.push_back(m_normals[values[randVal]]);
		values[randVal] = values[currentSize-1];
		currentSize--;
	}

	//Upsample
#pragma omp parallel for
	for (int i=0; i<(int)size() / int(divisionFactor); i++)
	{
		m_positions[i * divisionFactor] = newPositions[i];
		for (unsigned int j=1; j<divisionFactor; j++)
		{
			glm::vec3 randomVector(rand() / (float)RAND_MAX - 0.5, rand() / (float)RAND_MAX - 0.5, rand() / (float)RAND_MAX - 0.5);
			randomVector = glm::normalize(randomVector) * (rand() / (float)RAND_MAX)*0.1f;
			m_positions[i * divisionFactor + j] = (newPositions[i] + randomVector);
		}
	}
}

void Pn::clustering(vector<glm::vec3> &colors)
{
	colors.resize(size());
	//With Morton code
	unsigned power = 8;
	std::vector<std::pair<uint32_t, unsigned>> mortonCode(size());
#pragma omp parallel for
	for (int i=0; i<(int)size(); i++)
	{
		mortonCode[i].first = 0;
		mortonCode[i].second = i;
		uint32_t mortonValue1 = 0, mortonValue2 = 0, mortonValue3 = 0;
		mortonValue1 = (m_positions[i].x+2.f) / 4.f * pow(2, power);
		mortonValue2 = (m_positions[i].y+2.f) / 4.f * pow(2, power);
		mortonValue3 = (m_positions[i].z+2.f) / 4.f * pow(2, power);

		mortonCode[i].first = _pdep_u32(mortonValue1, 0x92492492) | _pdep_u32(mortonValue2, 0x49249249) | _pdep_u32(mortonValue3, 0x24924924);
	}
	std::sort(mortonCode.begin(), mortonCode.end(), [](std::pair<uint32_t, unsigned> val1, std::pair<uint32_t, unsigned> val2) {
		return val1.first < val2.first;
	});
	unsigned int counter = 0;
	uint32_t previousValue = mortonCode[0].first;
	glm::vec3 currentColor(rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
	vector<unsigned int> nbPointsPerCluster(1);
	nbPointsPerCluster[0] = 0;
	vector<unsigned int> clusters(size());
	for (unsigned int i=0; i<size(); i++)
	{
		if (previousValue != mortonCode[i].first)
		{
			counter++;
			nbPointsPerCluster.push_back(0);
			currentColor = glm::vec3(rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
		}
		nbPointsPerCluster[counter]++;
		previousValue = mortonCode[i].first;
//        colors[mortonCode[i].second] = currentColor;
		colors[i] = currentColor;
		clusters[i] = counter;
	}
	cout << counter+1 << " clusters" << endl;
}

void Pn::computePointsDensity()
{
	m_density.resize(size());
	//With Morton code
	unsigned power = pow(2, 5);
	std::vector<std::pair<uint32_t, unsigned>> mortonCode(size());
#pragma omp parallel for
	for (int i=0; i<(int)size(); i++)
	{
		mortonCode[i].first = 0;
		mortonCode[i].second = i;
		uint32_t mortonValue1 = 0, mortonValue2 = 0, mortonValue3 = 0;
		mortonValue1 = m_positions[i].x * power + pow(2,9);
		mortonValue2 = m_positions[i].y * power + pow(2,9);
		mortonValue3 = m_positions[i].z * power + pow(2,9);

		mortonCode[i].first = _pdep_u32(mortonValue1, 0x92492492) | _pdep_u32(mortonValue2, 0x49249249) | _pdep_u32(mortonValue3, 0x24924924);
	}
	std::sort(mortonCode.begin(), mortonCode.end(), [](std::pair<uint32_t, unsigned> val1, std::pair<uint32_t, unsigned> val2) {
		return val1.first < val2.first;
	});
	unsigned int counter = 0;
	uint32_t previousValue = mortonCode[0].first;
	glm::vec3 currentColor(rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
	vector<unsigned> nbPointsPerCluster(1);
	nbPointsPerCluster[0] = 0;
	vector<unsigned int> clusters(size());
	vector<vector<unsigned int> > ptsInClusters;
	vector<unsigned int> ptsInCurrentCluster;
	for (unsigned int i=0; i<size(); i++)
	{
		if (previousValue != mortonCode[i].first)
		{
			counter++;
			nbPointsPerCluster.push_back(0);
			ptsInClusters.push_back(ptsInCurrentCluster);
			ptsInCurrentCluster.clear();
		}
		nbPointsPerCluster[counter]++;
		ptsInCurrentCluster.push_back(mortonCode[i].second);
		previousValue = mortonCode[i].first;
		clusters[i] = counter;
	}
	ptsInClusters.push_back(ptsInCurrentCluster);
	cout << counter+1 << " clusters" << endl;

	vector<float> densityPerCluster(counter);

#pragma omp parallel for
	for (int i=0; i<(int)counter; i++)
	{
		float area = 0;
		if (ptsInClusters[i].size() > 2)
		{
			//Compute plane
			glm::vec3 n = glm::normalize(glm::cross((m_positions[ptsInClusters[i][1]] - m_positions[ptsInClusters[i][0]]), (m_positions[ptsInClusters[i][2]] - m_positions[ptsInClusters[i][0]])));
			float dl =glm::dot(n, m_positions[ptsInClusters[i][0]]);

			//compute cube
			glm::vec3 ptA = glm::vec3(floor(m_positions[ptsInClusters[i][0]].x * power) / (float)power,
					  floor(m_positions[ptsInClusters[i][0]].y * power) / (float)power,
					  floor(m_positions[ptsInClusters[i][0]].z * power) / (float)power );
			vector<glm::vec3> cubeVertices(8);
			cubeVertices[0] = ptA;
			cubeVertices[1] = glm::vec3(ptA.x + 1.f / (float)power, ptA.y, ptA.z);
			cubeVertices[2] = glm::vec3(ptA.x, ptA.y + 1.f / (float)power, ptA.z);
			cubeVertices[3] = glm::vec3(ptA.x, ptA.y, ptA.z + 1.f / (float)power);
			cubeVertices[4] = glm::vec3(ptA.x + 1.f / (float)power, ptA.y, ptA.z + 1.f / (float)power);
			cubeVertices[5] = glm::vec3(ptA.x + 1.f / (float)power, ptA.y + 1.f / (float)power, ptA.z);
			cubeVertices[6] = glm::vec3(ptA.x, ptA.y + 1.f / (float)power, ptA.z + 1.f / (float)power);
			cubeVertices[7] = glm::vec3(ptA.x + 1.f / (float)power, ptA.y + 1.f / (float)power, ptA.z + 1.f / (float)power);

			//compute polygon intersection
			vector<glm::vec3> polygonOfIntersection;
			glm::vec3 intersection;
			if (intersectEdge(dl, n, cubeVertices[0], cubeVertices[1], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[1], cubeVertices[4], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[4], cubeVertices[7], intersection))
				polygonOfIntersection.push_back(intersection);
			if (intersectEdge(dl, n, cubeVertices[1], cubeVertices[5], intersection))
				polygonOfIntersection.push_back(intersection);
			if (intersectEdge(dl, n, cubeVertices[0], cubeVertices[2], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[2], cubeVertices[5], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[5], cubeVertices[7], intersection))
				polygonOfIntersection.push_back(intersection);
			if (intersectEdge(dl, n, cubeVertices[2], cubeVertices[6], intersection))
				polygonOfIntersection.push_back(intersection);
			if (intersectEdge(dl, n, cubeVertices[0], cubeVertices[3], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[3], cubeVertices[6], intersection))
				polygonOfIntersection.push_back(intersection);
			else if (intersectEdge(dl, n, cubeVertices[6], cubeVertices[7], intersection))
				polygonOfIntersection.push_back(intersection);
			if (intersectEdge(dl, n, cubeVertices[3], cubeVertices[4], intersection))
				polygonOfIntersection.push_back(intersection);

			//compute intersection area
			unsigned nc = polygonOfIntersection.size();
			if (nc > 2)
			{
			glm::vec3 compute = glm::vec3(0, 0, 0);
			for (unsigned int j=0; j<nc; j++)
				compute += glm::cross(polygonOfIntersection[j], polygonOfIntersection[(j+1)%nc]);
			area = glm::dot(n, compute) / 2.f;
			area = fabs(area);
			}
			if (area < 0.001f)
				area = 1.f / (float)(power * power);
		}
		else
		{
			area = 1.f / (float)(power * power);
		}
		if (area == 0)
			area = 1.f / (float)(power * power);
		m_density[i] = area;
		densityPerCluster[i] = nbPointsPerCluster[i] / area;
	}
}
