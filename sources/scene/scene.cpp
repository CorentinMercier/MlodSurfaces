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

#include "scene.h"
#include "../interface/widgetopengl.h"
#include "../dual-contouring/dual_contouring.h"
#include "../dual-contouring/dual_contouringV2.h"
#include "upsample.h"


#ifdef COUT_DETAILED_STATISTICS
#include "../misc/DrawHistogram.h"
#endif



scene::scene()
{

}

scene::~scene()
{
	for (unsigned i=0; i<m_pnProjectedPoints.size(); i++)
		delete m_pnProjectedPoints[i];
}

void scene::reset()
{
	m_apss->eraseAndRestart();
}

void scene::compileShaders()
{
	ShaderInfo  shaders[] =
	{
		{ GL_VERTEX_SHADER, "../sources/opengl/shaders/shaderPoints.vert", 0 },
		{ GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shader.frag", 0 },
		{ GL_NONE, NULL, 0 }
	};
	m_basicProgram=loadShaders(shaders);
	if (m_basicProgram==0)
		cerr << "Error due to basic shaders" << endl;

	ShaderInfo  triangleShaders[] =
	{
		{ GL_VERTEX_SHADER, "../sources/opengl/shaders/shader.vert", 0 },
		{ GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shader.frag", 0 },
		{ GL_NONE, NULL, 0 }
	};
	m_triangleProgram=loadShaders(triangleShaders);
	if (m_triangleProgram==0)
		cerr << "Error due to triangle shaders" << endl;

	ShaderInfo singleColorShaders[] =
	{
		{ GL_VERTEX_SHADER, "../sources/opengl/shaders/shaderSingleColor.vert", 0 },
		{ GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shaderSingleColor.frag", 0 },
		{ GL_NONE, NULL, 0 }
	};
	m_singleColorProgram=loadShaders(singleColorShaders);
	if (m_singleColorProgram==0)
		cerr << "Error due to single color shaders" << endl;

	ShaderInfo perPointColorShaders[] =
	{
		{ GL_VERTEX_SHADER, "../sources/opengl/shaders/shaderWithPerPointColor.vert", 0 },
		{ GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shaderWithPerPointColor.frag", 0 },
		{ GL_NONE, NULL, 0 }
	};
	m_perPointColorProgram=loadShaders(perPointColorShaders);
	if (m_perPointColorProgram==0)
		cerr << "Error due to per point color shaders" << endl;

	//Lights definition
	glm::vec3 light1_position(22.0f, 16.0f, 50.0f);
	glm::vec3 light1_color(1.0f, 1.0f, 1.0f);
	float light1_intensity = 4.5f;
	glm::vec3 light2_position(-22.0f, 20.0f, 50.0f);
	glm::vec3 light2_color(1.0f, 0.8f, 1.0f);
	float light2_intensity = 4.5f;
	glm::vec3 ambient_color(0.3f, 0.3f, 0.3f);

	glUseProgram(m_basicProgram);
	glUniform3f (glGetUniformLocation (m_basicProgram, "lightSource1.position"), light1_position.x, light1_position.y, light1_position.z);
	glUniform3f (glGetUniformLocation (m_basicProgram, "lightSource1.color"), light1_color[0], light1_color[1], light1_color[2]);
	glUniform1f (glGetUniformLocation (m_basicProgram, "lightSource1.intensity"), light1_intensity);
	glUniform3f (glGetUniformLocation (m_basicProgram, "lightSource2.position"), light2_position.x, light2_position.y, light2_position.z);
	glUniform3f (glGetUniformLocation (m_basicProgram, "lightSource2.color"), light2_color[0], light2_color[1], light2_color[2]);
	glUniform1f (glGetUniformLocation (m_basicProgram, "lightSource2.intensity"), light2_intensity);
	glUniform3f (glGetUniformLocation (m_basicProgram, "ambient_color"), ambient_color[0], ambient_color[1], ambient_color[2]);

	glUseProgram(m_triangleProgram);
	glUniform3f (glGetUniformLocation (m_triangleProgram, "lightSource1.position"), light1_position.x, light1_position.y, light1_position.z);
	glUniform3f (glGetUniformLocation (m_triangleProgram, "lightSource1.color"), light1_color[0], light1_color[1], light1_color[2]);
	glUniform1f (glGetUniformLocation (m_triangleProgram, "lightSource1.intensity"), light1_intensity);
	glUniform3f (glGetUniformLocation (m_triangleProgram, "lightSource2.position"), light2_position.x, light2_position.y, light2_position.z);
	glUniform3f (glGetUniformLocation (m_triangleProgram, "lightSource2.color"), light2_color[0], light2_color[1], light2_color[2]);
	glUniform1f (glGetUniformLocation (m_triangleProgram, "lightSource2.intensity"), light2_intensity);
	glUniform3f (glGetUniformLocation (m_triangleProgram, "ambient_color"), ambient_color[0], ambient_color[1], ambient_color[2]);

	glUseProgram(m_perPointColorProgram);
	glUniform3f (glGetUniformLocation (m_perPointColorProgram, "lightSource.position"), light1_position.x, light1_position.y, light1_position.z);
	glUniform3f (glGetUniformLocation (m_perPointColorProgram, "lightSource.color"), light1_color[0], light1_color[1], light1_color[2]);
	glUniform1f (glGetUniformLocation (m_perPointColorProgram, "lightSource.intensity"), light1_intensity);
	glUniform3f (glGetUniformLocation (m_perPointColorProgram, "ambient_color"), ambient_color[0], ambient_color[1], ambient_color[2]);
}

void scene::createBuffers()
{
	glGenVertexArrays(1, &m_vao); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboPointSet); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vao2); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboPointSetToProject); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vao3); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboPointSetProjected); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vao5); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboWithColors); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vaoSelection); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboSelection); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vaoDualContouring); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboDualContouring); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_iboDualContouring); PRINT_OPENGL_ERROR();
	glGenVertexArrays(1, &m_vaoComparison); PRINT_OPENGL_ERROR();
	glCreateBuffers(1, &m_vboComparison); PRINT_OPENGL_ERROR();
}

void scene::eraseBuffers()
{
	glDeleteVertexArrays(1, &m_vao); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboPointSet); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vao2); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboPointSetToProject); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vao3); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboPointSetProjected); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vao5); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboWithColors); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vaoSelection); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboSelection); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vaoDualContouring); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboDualContouring); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_iboDualContouring); PRINT_OPENGL_ERROR();
	glDeleteVertexArrays(1, &m_vaoComparison); PRINT_OPENGL_ERROR();
	glDeleteBuffers(1, &m_vboComparison); PRINT_OPENGL_ERROR();
}

void scene::generate_adaptive_octree()
{
	OctreeGenerator octree_gen(m_apss, *m_weightKernel);
	m_totalUpsamplingNodes.clear();
	cout << "Computing adaptive octree for upsampling ..." << endl;
	CPUTimer t;
	t.start();
	m_adaptiveOctree = octree_gen.generate(m_maxDepth, &m_totalUpsamplingNodes/*, transform*/);
	t.stop();
	cout << "   [OK] " << m_totalUpsamplingNodes.size() << " nodes were crossed at depth " << m_maxDepth << endl;
	m_nbOfCrossedNodes = m_totalUpsamplingNodes.size();
	m_adaptiveOctreeTime = t.printElapsed("Time for adaptive octree creation");
}

void scene::generate_usual_octree()
{
	OctreeGenerator octree_gen(m_apss, *m_weightKernel);
	cout << "Computing usual octree for upsampling ..." << endl;
	CPUTimer t;
	t.start();
	m_usualOctree = octree_gen.generateFromPoints(m_maxDepth, m_originalPn);
	cout << "   [OK] " << endl;
	m_usualOctreeTime = t.printElapsed("Time for the creation of the octree on the input pointset");
}

size_t scene::count_nodes_on_projected_points()
{
	OctreeGenerator octree_gen(m_apss, *m_weightKernel);
	std::vector<octree_id> nodes;

	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	const Eigen::AlignedBox3f B = octzone.bounds();
	Eigen::Matrix4f unit2bbox = Eigen::Matrix4f::Zero();
	unit2bbox(0, 0) = B.diagonal()[0];
	unit2bbox(1, 1) = B.diagonal()[1];
	unit2bbox(2, 2) = B.diagonal()[2];
	unit2bbox(0, 3) = B.min()[0];
	unit2bbox(1, 3) = B.min()[1];
	unit2bbox(2, 3) = B.min()[2];
	unit2bbox(3, 3) = 1;

	cout << "Uniformizing tree ..." << endl;

	sorted_octree<char> totalUpsamplingOctree;
	glm::vec3 new_diag(1), new_origin(0);
	CPUTimer t;
	t.start();
	if (m_useAdaptiveOctree)
		uniformize(m_maxDepth, m_totalUpsamplingNodes, totalUpsamplingOctree, new_diag, new_origin, m_cut);
	else
		uniformizeFromPoints(m_usualOctree, totalUpsamplingOctree, new_diag, new_origin);
	m_uniformizeTime = t.printElapsedAndReset("[ OK ] ");
	cout << totalUpsamplingOctree.size() << " nodes" << endl;

	Eigen::Matrix4f sortedOctreeToUnit = Eigen::Matrix4f::Zero();
	sortedOctreeToUnit(0, 0) = new_diag[0];
	sortedOctreeToUnit(1, 1) = new_diag[1];
	sortedOctreeToUnit(2, 2) = new_diag[2];
	sortedOctreeToUnit(0, 3) = new_origin[0];
	sortedOctreeToUnit(1, 3) = new_origin[1];
	sortedOctreeToUnit(2, 3) = new_origin[2];
	sortedOctreeToUnit(3, 3) = 1;

	Eigen::Matrix4f new2bbox = unit2bbox * sortedOctreeToUnit;

	Eigen::Matrix4f bbox2orig = Eigen::Matrix4f::Zero();
	bbox2orig(0, 0) = m_originalPn.getNorme();
	bbox2orig(1, 1) = m_originalPn.getNorme();
	bbox2orig(2, 2) = m_originalPn.getNorme();
	bbox2orig(0, 3) = m_originalPn.getOffset().x;
	bbox2orig(1, 3) = m_originalPn.getOffset().y;
	bbox2orig(2, 3) = m_originalPn.getOffset().z;
	bbox2orig(3, 3) = 1;
	const Eigen::Matrix4f transform = bbox2orig * unit2bbox * sortedOctreeToUnit;

	Eigen::Vector4f min = new2bbox * Eigen::Vector4f(  0.0f,  0.0f,  0.0f, 1.f);
	Eigen::Vector4f max = new2bbox * Eigen::Vector4f(  1.0f,  1.0f,  1.0f, 1.f);

	nodes = octree_gen.generateFromPoints(m_maxDepth, *m_pnProjectedPoints[0], Eigen::AlignedBox3f(Eigen::Vector3f(min.x(), min.y(), min.z()), Eigen::Vector3f(max.x(), max.y(), max.z())));
	regularize_up(nodes);
	save_octree_ply("Octree.ply", nodes, node_color_by_depth, transform);
	save_octree_ply("octree.ply", m_totalUpsamplingNodes, node_color_by_depth, bbox2orig * unit2bbox);
	cout << "Nodes on projected points: " << nodes.size() << endl;
	return nodes.size();
}

void scene::initBuffers()
{
	cout << "Buffer initialization ..." << endl;
	m_selection.resize(0);

	//Fill GPU memory with octree
	m_apss->copyToGPU(*m_weightKernel, m_nbOfPoints);

	generate_adaptive_octree();
	generate_usual_octree();

	// Create points based on the octree
	Pn newPn(generateFromAdaptiveOctree());
	m_pnToProject2 = newPn;
	m_pnToProject2.order();
	// Create random points to project
	glm::vec3 mini = m_pn.getMini() + glm::vec3(m_offset[0], m_offset[1], m_offset[2]);
	glm::vec3 maxi = m_pn.getMaxi() + glm::vec3(m_offset[3], m_offset[4], m_offset[5]);
	m_pnToProject.initialize(m_nbOfPoints, mini, maxi, m_sizeOfGrid);

	m_apss->copyPointsToGPU(m_pnToProject2.size(), m_pnToProject2.positionData());


	//*****************************************//
	// Create OPENGL VBO for the pointsets      //
	//*****************************************//
	glEnable(GL_POINT_SMOOTH); PRINT_OPENGL_ERROR();
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE); PRINT_OPENGL_ERROR();
	glEnable(GL_POLYGON_OFFSET_LINE);//For wireframe mode

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSet);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size()*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_pn.size(), m_pn.positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size(), 3*sizeof(float)*m_pn.size(), m_pn.normalData());PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSet);
	GLint loc1=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc1);
	GLint loc2=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pn.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc2);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetToProject);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size(), m_pnToProject.positionData(), GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao2);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetToProject);
	GLint loc3=glGetAttribLocation(m_singleColorProgram, "position");
	glVertexAttribPointer(loc3, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc3);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size()*2*m_maxUpsamples, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao3);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);
	GLint loc4=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc4, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc4);
	GLint loc5=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc5, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size() * m_maxUpsamples));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc5);


	glBindBuffer(GL_ARRAY_BUFFER, m_vboWithColors);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size()*3, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao5);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboWithColors);
	GLint loc7=glGetAttribLocation(m_perPointColorProgram, "position");
	glVertexAttribPointer(loc7, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc7);
	GLint loc8=glGetAttribLocation(m_perPointColorProgram, "normal");
	glVertexAttribPointer(loc8, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size())); PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc8);
	GLint loc9=glGetAttribLocation(m_perPointColorProgram, "perPointColor");
	glVertexAttribPointer(loc9, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size() * 2));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc9);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size()*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vaoSelection);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	GLint loc10=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc10, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc10);
	GLint loc11=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc11, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pn.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc11);

	initUpsampling();

	cout << "Buffer initialization finished" << endl << endl;
}

double scene::updateBuffers(bool orderWithMorton)
{
	m_pnToProject2.initialize(m_nbOfPoints, m_pn, m_pwidget->getRadiusFactor(), m_sizeOfGrid, orderWithMorton);
	glm::vec3 mini = m_pn.getMini() + glm::vec3(m_offset[0], m_offset[1], m_offset[2]);
	glm::vec3 maxi = m_pn.getMaxi() + glm::vec3(m_offset[3], m_offset[4], m_offset[5]);
	m_pnToProject.initialize(m_nbOfPoints, mini, maxi, m_sizeOfGrid, orderWithMorton);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSet);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size()*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_pn.size(), m_pn.positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size(), 3*sizeof(float)*m_pn.size(), m_pn.normalData());PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSet);
	GLint loc1=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc1);
	GLint loc2=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pn.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc2);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetToProject);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size(), m_pnToProject.positionData(), GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao2);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetToProject);
	GLint loc3=glGetAttribLocation(m_singleColorProgram, "position");
	glVertexAttribPointer(loc3, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc3);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size()*2*m_maxUpsamples, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao3);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);
	GLint loc4=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc4, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc4);
	GLint loc5=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc5, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size() * m_maxUpsamples));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc5);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboWithColors);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnToProject.size()*3, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vao5);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboWithColors);
	GLint loc7=glGetAttribLocation(m_perPointColorProgram, "position");
	glVertexAttribPointer(loc7, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc7);
	GLint loc8=glGetAttribLocation(m_perPointColorProgram, "normal");
	glVertexAttribPointer(loc8, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size())); PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc8);
	GLint loc9=glGetAttribLocation(m_perPointColorProgram, "perPointColor");
	glVertexAttribPointer(loc9, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnToProject.size() * 2));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc9);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size()*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vaoSelection);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	GLint loc10=glGetAttribLocation(m_basicProgram, "position");
	glVertexAttribPointer(loc10, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc10);
	GLint loc11=glGetAttribLocation(m_basicProgram, "normal");
	glVertexAttribPointer(loc11, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pn.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc11);

	double mortonTime = resetPoints(true);

	cout << "Buffers updated" << endl;
	return mortonTime;
}

bool scene::load_scene(string filename)
{
	compileShaders();

	//*****************************************//
	// Load pointset                           //
	//*****************************************//

	auto startLoading = chrono::steady_clock::now();
	auto start = chrono::steady_clock::now();

	// Load pointset data
	cout << "  Load data from file " << filename << " ... " << endl;
	m_pn.load(filename, m_pwidget->getRadiusFactor());
	cout << m_pn.size() << " points found -- [OK] " << endl;
	m_originalPn = m_pn;

	// Create octree
	m_apssOctree = new APSSOctree();
	m_apssOctree->build(m_pn.getPositions(), m_pn.getNormals());
	cout << "Nb of leaves in root: " << m_apssOctree->get_root()->nbOfNodes() << endl;

	OctreeGenerator octree_gen(m_apss, *m_weightKernel);
	m_usualOctree = octree_gen.generateFromPoints(m_maxDepth, m_pn);

#ifdef COUT_DETAILED_STATISTICS
	debug_kdtree.build( m_pn.getPositions() );
	debug_kdtree.index_nodes_by_breadth_first_search();
#endif

	//Initialize weight function
	m_weightKernel = new Kernel();
	if(false) {
		m_weightKernel->sigma = 0.01 * m_apssOctree->getBoundingBox().radius();
		m_weightKernel->type = Kernel::SINGULAR;
		m_weightKernel->exponent = 3.0;
		m_weightKernel->constant_shift = 0.0000000000001;
	}
	else {
		m_weightKernel->type = Kernel::GAUSSIAN_MULTIPLE;
		m_weightKernel->sigmas.clear();
		m_weightKernel->sigmas.push_back(0.02 * m_apssOctree->getBoundingBox().radius());
		m_weightKernel->sigmas.push_back(0.05 * m_apssOctree->getBoundingBox().radius());
		m_weightKernel->sigmas.push_back(0.1 * m_apssOctree->getBoundingBox().radius());
		m_weightKernel->sigmas.push_back(0.3 * m_apssOctree->getBoundingBox().radius());
		m_weightKernel->sigmas.push_back(1.5 * m_apssOctree->getBoundingBox().radius());
		m_weightKernel->constant_shift = 0.0;
	}

	//Load parameters if any exist
	if (!loadParameters(filename))
		m_apss = new Fastapss(m_apssOctree);
	m_pwidget->connectKernel();

	auto end = chrono::steady_clock::now();
	auto diff = end-start;
	cout << "Tree time : " << chrono::duration <double, milli> (diff).count() << " ms" << endl;

	createBuffers();
	cout << "Filling buffers ..." << endl;
	initBuffers();
	cout << "   [OK] " << endl;

	cout << "Run drawing" << endl;

	auto endLoading = chrono::steady_clock::now();
	auto diffLoading = endLoading-startLoading;
	cout << "Total loading time : " << chrono::duration <double, milli> (diffLoading).count() << " ms" << endl;
	m_time.start();
	return true;
}

void scene::changeMethod(bool fast)
{
	if (m_method == fast) return;
	float cutOff = m_apss->getScalingProtectionSphere();
	delete m_apss;
	if (fast)
		m_apss = new Fastapss(m_apssOctree, *m_weightKernel, m_nbOfPoints);
	else
		m_apss = new GlobalAPSS(m_apssOctree, *m_weightKernel, m_nbOfPoints);
	m_apss->setMinDepth(m_minDepth);
	m_apss->setScalingProtectionSphere(cutOff);
	if (m_pnProjectedPoints[m_currentProjection]->size() == 0)
	{
		delete m_pnProjectedPoints[m_currentProjection];
		if (m_random)
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject);
		else
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject2);
	}
	m_apss->copyPointsToGPU(m_nbOfPoints, m_pnProjectedPoints[m_currentProjection]->positionData());
	m_method = fast;
}

void scene::save_octree(const std::string& filename, unsigned int min_depth, unsigned int max_depth, bool use_empty_triangles) const
{
	size_t num_lines = 0;
	std::function<void(const octreeNode*)> count_lines = [&](const octreeNode* node)
	{
		if(node->depth >= min_depth)
			num_lines += 12;
		if(node->depth < max_depth)
			for(unsigned int c = 0; c < 8; c++)
				if(node->children[c]) count_lines(node->children[c]);
	};
	count_lines(m_apssOctree->get_root());
	const size_t num_vertices = num_lines * (use_empty_triangles ? 3 : 2);

	std::ofstream out(filename);
	out << "ply\n"
		<< "format ascii 1.0\n"
		<< "element vertex " << num_vertices << "\n"
		<< "property float x\n"
		<< "property float y\n"
		<< "property float z\n"
		<< "property uchar red\n"
		<< "property uchar green\n"
		<< "property uchar blue\n";
	if(use_empty_triangles)
		out << "element face " << num_lines << "\n"
		<< "property list uchar int vertex_index\n";
	else
		out << "element edge " << num_lines << "\n"
		<< "property int vertex1\n"
		<< "property int vertex2\n";
	out << "end_header\n";

	auto node_color_by_depth = [&](unsigned int depth)
	{
		const glm::ivec3 cols[] =
		{
			{ 0,   255, 0   },
			{ 128, 255, 128 },
			{ 128, 255, 255 },
			{ 0,   0,   255 },
			{ 128, 128, 255 },
			{ 255, 128, 255 },
			{ 255, 0,   0   },
			{ 255, 128, 128 },
			{ 255, 255, 128 }
		};
		const size_t num_colors = (sizeof(cols) / sizeof(cols[0]));
		return cols[depth % num_colors];
	};

	auto add_vert = [&](const glm::vec3& pos, const glm::ivec3& col)
	{
		out << pos[0] << ' ' << pos[1] << ' ' << pos[2] << ' '
			<< col[0] << ' ' << col[1] << ' ' << col[2] << '\n';
	};
	auto add_line = [&](const glm::ivec3& color, const glm::vec3& a, const glm::vec3& b)
	{
		add_vert(a, color);
		add_vert(b, color);
		if(use_empty_triangles) add_vert((a + b) / 2.0f, color);
	};
	std::function<void(const octreeNode*)> traverse = [&](const octreeNode* node)
	{
		if(node->depth >= min_depth)
		{
			const glm::ivec3 color = node_color_by_depth(node->depth);
			add_line(color, node->boundingBox.corner(0), node->boundingBox.corner(1));
			add_line(color, node->boundingBox.corner(0), node->boundingBox.corner(2));
			add_line(color, node->boundingBox.corner(2), node->boundingBox.corner(3));
			add_line(color, node->boundingBox.corner(1), node->boundingBox.corner(3));
			add_line(color, node->boundingBox.corner(4), node->boundingBox.corner(5));
			add_line(color, node->boundingBox.corner(4), node->boundingBox.corner(6));
			add_line(color, node->boundingBox.corner(6), node->boundingBox.corner(7));
			add_line(color, node->boundingBox.corner(5), node->boundingBox.corner(7));
			add_line(color, node->boundingBox.corner(0), node->boundingBox.corner(4));
			add_line(color, node->boundingBox.corner(1), node->boundingBox.corner(5));
			add_line(color, node->boundingBox.corner(2), node->boundingBox.corner(6));
			add_line(color, node->boundingBox.corner(3), node->boundingBox.corner(7));
		}
		if(node->depth < max_depth)
			for(unsigned int c = 0; c < 8; c++)
				if(node->children[c]) traverse(node->children[c]);
	};
	traverse(m_apssOctree->get_root());

	for(size_t l = 0; l < num_lines; l++)
	{
		const size_t c = use_empty_triangles ? 3 : 2;
		if(c != 2) out << c << ' ';
		for(size_t v = 0; v < c; v++)
		{
			if(v != 0) out << ' ';
			out << (l * c + v);
		}
		out << '\n';
	}
	std::cout << "Saved octree (depths: min = " << min_depth << ", max = " << max_depth << ") to '" << filename << "'\n";
}

void scene::draw_scene()
{
	glUseProgram(m_basicProgram);
	glUniformMatrix4fv(get_uni_loc(m_basicProgram,"camera_modelview"),1,false, glm::value_ptr(m_pwidget->cam.getModelview()));  PRINT_OPENGL_ERROR();
	glUniformMatrix4fv(get_uni_loc(m_basicProgram,"camera_projection"),1,false, glm::value_ptr(m_pwidget->cam.getProjection()));   PRINT_OPENGL_ERROR();
	glUniformMatrix4fv(get_uni_loc(m_basicProgram,"normal_matrix"),1,false, glm::value_ptr(m_pwidget->cam.getNormal()));   PRINT_OPENGL_ERROR();
	glUniform1i(get_uni_loc(m_basicProgram,"normal_mode"), normalMode);   PRINT_OPENGL_ERROR();

	glUniform1i(get_uni_loc(m_basicProgram, "pointSize"), m_ptSize);

	if (renderSelection && m_selection.size() != 0)
	{
		glUniform4f(get_uni_loc(m_basicProgram, "pointsColor"), 0.8f, 0.f, 0.2f, 1.f);
		glBindVertexArray(m_vaoSelection);
		glDrawArrays(GL_POINTS, 0, m_selection.size());
	}

	if (renderPointSet)
	{
		glUniform4f(get_uni_loc(m_basicProgram, "pointsColor"), 1.f, 1.f, 1.f, 1.f);
		glBindVertexArray(m_vao);
		glDrawArrays(GL_POINTS, 0, m_pn.size());
	}

	if (renderKernel)
	{
		glUniform1i(get_uni_loc(m_basicProgram, "pointSize"), 10);
		glUniform4f(get_uni_loc(m_basicProgram, "pointsColor"), 0.2941f, 0.f, 0.5098f, 1.f);
		m_weightKernel->draw(m_basicProgram);
		glUniform1i(get_uni_loc(m_basicProgram, "pointSize"), m_ptSize);
	}

	if (renderProjectedPoints && m_pnProjectedPoints[0]->size() != 0)
	{
		glUniform4f(get_uni_loc(m_basicProgram, "pointsColor"), 0.f, 1.f, 0.05f, 1.f);
		glBindVertexArray(m_vao3);
		if (m_nbOfInvalidPts == 0)
			glDrawArrays(GL_POINTS, 0, m_pnProjectedPoints[0]->size() * (min(m_currentProjection + 1, m_maxUpsamples)));
		else
		{
			for (unsigned i=0; i<(min(m_currentProjection + 1, m_maxUpsamples)); i++)
				glDrawArrays(GL_POINTS, m_nbOfPoints * i, m_pnProjectedPoints[i]->validSize());
		}
	}

	if (renderDualContour && m_nbOfIndices != 0)
	{
		glUseProgram(m_triangleProgram);
		glUniformMatrix4fv(get_uni_loc(m_triangleProgram,"camera_modelview"),1,false, glm::value_ptr(m_pwidget->cam.getModelview()));  PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_triangleProgram,"camera_projection"),1,false, glm::value_ptr(m_pwidget->cam.getProjection()));   PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_triangleProgram,"normal_matrix"),1,false, glm::value_ptr(m_pwidget->cam.getNormal()));   PRINT_OPENGL_ERROR();
		glUniform1i(get_uni_loc(m_triangleProgram,"normal_mode"), normalMode);   PRINT_OPENGL_ERROR();
		glUniform4f(get_uni_loc(m_triangleProgram,"pointsColor"), 0.8f, 0.1f, 0.1f, 1.f);   PRINT_OPENGL_ERROR();
		glBindVertexArray(m_vaoDualContouring);
		glDrawElements(GL_TRIANGLES, m_nbOfIndices, GL_UNSIGNED_INT, (void*)0); PRINT_OPENGL_ERROR();
		if (wireframe)
		{
			glEnable(GL_POLYGON_OFFSET_LINE);//For wireframe mode
			glPolygonOffset(-1.0, -1.0);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glUniform4f(get_uni_loc(m_triangleProgram,"pointsColor"), 0.0f, 0.0f, 0.0f, 1.f);   PRINT_OPENGL_ERROR();
			glDrawElements(GL_TRIANGLES, m_nbOfIndices, GL_UNSIGNED_INT, (void*)0); PRINT_OPENGL_ERROR();
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
	}

	glUseProgram(m_singleColorProgram);
	glUniformMatrix4fv(get_uni_loc(m_singleColorProgram,"camera_modelview"),1,false, glm::value_ptr(m_pwidget->cam.getModelview()));  PRINT_OPENGL_ERROR();
	glUniformMatrix4fv(get_uni_loc(m_singleColorProgram,"camera_projection"),1,false, glm::value_ptr(m_pwidget->cam.getProjection()));   PRINT_OPENGL_ERROR();
	glUniform1i(get_uni_loc(m_singleColorProgram, "pointSize"), m_ptSize);

	if (renderPointsToProject)
	{
		glUniform3f(get_uni_loc(m_singleColorProgram,"pointsColor"), 0.8f, 0.1f, 0.1f);   PRINT_OPENGL_ERROR();
		glBindVertexArray(m_vao2);
		glDrawArrays(GL_POINTS, 0, m_pnToProject.size());
	}

	if (renderInvalidPoints && m_nbOfInvalidPts != 0)
	{
		glUniform3f(get_uni_loc(m_singleColorProgram, "pointsColor"), 0.9f, 0.2f, 0.15f); PRINT_OPENGL_ERROR();
		glBindVertexArray(m_vao3);
		for (unsigned i=0; i<(min(m_currentProjection + 1, m_maxUpsamples)); i++)
			if (m_pnProjectedPoints[i]->size() - m_pnProjectedPoints[i]->validSize() > 0)
				glDrawArrays(GL_POINTS, m_nbOfPoints * i + m_pnProjectedPoints[i]->validSize(), m_pnProjectedPoints[i]->size() - m_pnProjectedPoints[i]->validSize());
	}

	if (renderPerPointColor)
	{
		glUseProgram(m_perPointColorProgram);
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"camera_modelview"),1,false, glm::value_ptr(m_pwidget->cam.getModelview()));  PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"camera_projection"),1,false, glm::value_ptr(m_pwidget->cam.getProjection()));   PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"normal_matrix"),1,false, glm::value_ptr(m_pwidget->cam.getNormal()));   PRINT_OPENGL_ERROR();
		glUniform1i(get_uni_loc(m_perPointColorProgram, "pointSize"), m_ptSize);
		glBindVertexArray(m_vao5);
		glDrawArrays(GL_POINTS, 0, m_pnProjectedPoints[m_currentProjection]->size());
	}

	if (renderComparison && m_pnForComparison.size() != 0)
	{
		glUseProgram(m_perPointColorProgram);
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"camera_modelview"),1,false, glm::value_ptr(m_pwidget->cam.getModelview()));  PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"camera_projection"),1,false, glm::value_ptr(m_pwidget->cam.getProjection()));   PRINT_OPENGL_ERROR();
		glUniformMatrix4fv(get_uni_loc(m_perPointColorProgram,"normal_matrix"),1,false, glm::value_ptr(m_pwidget->cam.getNormal()));   PRINT_OPENGL_ERROR();
		glUniform1i(get_uni_loc(m_perPointColorProgram, "pointSize"), m_ptSize);
		glBindVertexArray(m_vaoComparison);
		glDrawArrays(GL_POINTS, 0, m_pnForComparison.size());
	}

	glFlush();
	m_pwidget->swapBuffers();

	m_time.stop();
	m_time.start();
}

void scene::set_widget(WidgetOpenGL *widget_param)
{
	m_pwidget=widget_param;
}

Pn scene::generateFromOctree(const std::vector<octree_id> &octree)
{
	std::vector<glm::vec3> newPoints(m_nbOfPoints);
	std::vector<glm::vec3> newNormals(m_nbOfPoints);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<unsigned int> dist(0, octree.size() - 1);
	std::uniform_real_distribution<float> dist2(-0.5f, 0.5f);
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	for (unsigned i=0; i<m_nbOfPoints; i++)
	{
		const unsigned int randomPoint = dist(mt);
		const octree_id & octId = octree[randomPoint];
		const float randomX = dist2(mt), randomY = dist2(mt), randomZ = dist2(mt);
		const float move = 4.0f / float(octree_id::width(octId.depth));// Times 2 because we are between -1, -1, -1 and 1, 1, 1, times 2 again for better smoothing, 4->becomes 2 for histograms on points in the bounding box
		Eigen::Vector3f c = octzone.center(octId);
		glm::vec3 center = glm::vec3(c.x(), c.y(), c.z());
		newPoints[i] = center + glm::vec3(randomX * move, randomY * move, randomZ * move);
		newNormals[i] = newPoints[i];
	}
	return Pn(std::move(newPoints), std::move(newNormals), m_originalPn.getOffset(), m_originalPn.getNorme());
}

Pn scene::generateFromAdaptiveOctree()
{
	return generateFromOctree(m_adaptiveOctree);
}

Pn scene::generateFromUsualOctree()
{
	return generateFromOctree(m_usualOctree);
}

Pn scene::generateFromBoundingBox()
{
	std::vector<octree_id> bbox;
	bbox.push_back(octree_id() );
	return generateFromOctree(bbox);
}

bool scene::upsample()
{
	double t;
	return upsample(t);
}

bool scene::upsample(double &mortonCodeTime)
{
	if (m_currentProjection + 1 < m_maxUpsamples)
	{
		Pn * newPnToProject;
		if (m_useAdaptiveOctree)
		{
			newPnToProject = new Pn(generateFromAdaptiveOctree());
			mortonCodeTime = newPnToProject->order();
		}
		else
		{
			newPnToProject = new Pn(generateFromUsualOctree());
			mortonCodeTime = newPnToProject->order();
		}
		m_currentProjection++;
		m_pnProjectedPoints.push_back(newPnToProject);
		m_apss->erasePointsFromGPU();
		m_apss->copyPointsToGPU(m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());
		return true;
	}
	else
	{
		m_currentProjection = m_maxUpsamples;
		return false;
	}
}


void scene::output_traversed_nodes_statistics() {
	unsigned int grid_size_here = 16;

#ifdef COUT_DETAILED_STATISTICS
	std::vector< double > histo_data;

	std::vector<unsigned int> knnToQuery = {20, 40, 100, 200, 500};
	std::vector<std::vector< double > > histo_data_knn(knnToQuery.size());
#endif

	Pn * newPn = m_pnProjectedPoints[m_currentProjection];
	std::vector< std::pair<uint32_t, unsigned> >  mortonCode;
	newPn->computeMortonOrder( grid_size_here , mortonCode );

	std::map< uint32_t , std::set< unsigned int > > nodesPerMortonCase;
	std::map< uint32_t , unsigned int > pointsPerMortonCase;
	std::map< uint32_t , unsigned int > nodesLazyPerMortonCase;

	for( auto c : mortonCode ) {
		uint32_t code = c.first;
		unsigned int pIdx = c.second;
		glm::vec3 q = (*newPn)[ pIdx ];
		std::set< unsigned int > traversed_nodes;
		m_apssOctree->get_root()->gatherTraversedNodes(q , traversed_nodes , m_minDepth , m_apssOctree->getScalingProtectionSphere() );

#ifdef COUT_DETAILED_STATISTICS
		histo_data.push_back( traversed_nodes.size() );

		for (unsigned i=0; i<knnToQuery.size(); i++)
		{
			kdTreeNN::searchResults knnRes;
			std::set< unsigned int > traversed_nodes_for_knn;
			debug_kdtree.gatherNodesForKnearest( q , knnToQuery[i] , traversed_nodes_for_knn , knnRes );
			histo_data_knn[i].push_back( traversed_nodes_for_knn.size() );
		}
#endif

		pointsPerMortonCase[ code ]++;
		nodesLazyPerMortonCase[ code ] += traversed_nodes.size();
		nodesPerMortonCase[ code ].insert( traversed_nodes.begin() , traversed_nodes.end() );
	}

#ifdef COUT_DETAILED_STATISTICS
	{
		unsigned int data_max_value = histo_data[0];
		unsigned int data_min_value = histo_data[0];
		for( unsigned int d : histo_data ) {
			data_max_value = std::max< unsigned int > ( data_max_value , d );
			data_min_value = std::min< unsigned int > ( data_min_value , d );
		}

		DrawHistogram histo;
		histo.setBackgroundColor();
		histo.fillData(histo_data , 512);
		histo.buildSceneFromData(512,512 , QColor( 120 , 120 , 200 ));
		histo.renderImage(QString("nb_nodes_ours_cutoff_%1-%2-%3.png").arg(this->m_apssOctree->getScalingProtectionSphere()).arg(data_min_value).arg(data_max_value).toStdString(),512,512);

		fstream file;
		file.open("histo.txt", ios_base::out);
		for (unsigned i=0; i<histo_data.size(); i++)
			file <<  histo_data[i] << endl;
		file.close();
	}
	{
		for (unsigned i=0; i<knnToQuery.size(); i++)
		{
			unsigned int data_max_value = histo_data_knn[i][0];
			unsigned int data_min_value = histo_data_knn[i][0];
			for( unsigned int d : histo_data_knn[i] ) {
				data_max_value = std::max< unsigned int > ( data_max_value , d );
				data_min_value = std::min< unsigned int > ( data_min_value , d );
			}

			DrawHistogram histo;
			histo.setBackgroundColor();
			histo.fillData(histo_data_knn[i] , 512);
			histo.buildSceneFromData(512,512 , QColor( 200 , 120 , 120 ));
			histo.renderImage(QString("nb_nodes_knn_cutoff_%1-%2-%3_%4.png").arg(this->m_apssOctree->getScalingProtectionSphere()).arg(data_min_value).arg(data_max_value).arg(knnToQuery[i]).toStdString(),512,512);

			fstream file;
			file.open("histoknn" + to_string(knnToQuery[i]) + ".txt", ios_base::out);
			for (unsigned j=0; j<histo_data_knn[i].size(); j++)
				file <<  histo_data_knn[i][j] << endl;
			file.close();
		}
	}
#endif

	for( std::map< uint32_t , std::set< unsigned int > >::const_iterator it = nodesPerMortonCase.begin() ; it != nodesPerMortonCase.end() ; ++it ) {
		std::cout << it->second.size() << " different nodes were traversed for " << pointsPerMortonCase[it->first]
				  << " projection points (grid size = " << grid_size_here << "^3), total with duplicates = " << nodesLazyPerMortonCase[it->first] << std::endl;
	}
}

void scene::initUpsampling()
{
	for (unsigned i=0; i<m_pnProjectedPoints.size(); i++)
		delete m_pnProjectedPoints[i];
	m_pnProjectedPoints.resize(1);
	m_pnProjectedPoints[0] = new Pn();
	m_currentProjection = 0;
	m_nbOfInvalidPts = 0;
}

double scene::projectPoints(bool updateBuffer)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
	//Create object for projected points
	if (m_pnProjectedPoints[m_currentProjection]->size() == 0)
	{
		delete m_pnProjectedPoints[m_currentProjection];
		if (m_random)
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject);
		else
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject2);
	}
	start = std::chrono::high_resolution_clock::now();
	m_apss->project(m_pnProjectedPoints[m_currentProjection]->validSize(), m_pnProjectedPoints[m_currentProjection]->positionData(), m_pnProjectedPoints[m_currentProjection]->normalData(), 1, *m_weightKernel);
	stop = std::chrono::high_resolution_clock::now();

#ifdef CHECK
		unsigned nbOfNan;
		checkForInvalidPoints(nbOfNan);
		if (nbOfNan != 0)
		{
			cerr << "--------------------------" << endl;
			cerr << "Warning, " << nbOfNan << " NaN values found" << endl;
			cerr << "This can happen when constant shift = 0" << endl;
			cerr << "Or if the first sigmas are too small" << endl;
			cerr << "--------------------------" << endl;
		}
#endif

	if (updateBuffer)
	{
		glBindVertexArray(m_vao3);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_currentProjection * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + m_currentProjection*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->normalData());PRINT_OPENGL_ERROR();
	}
	else if (m_currentProjection != 0) //Put last pointset in buffer
	{
		glBindVertexArray(m_vao3);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_currentProjection) * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection - 1]->size(), m_pnProjectedPoints[m_currentProjection - 1]->positionData());PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + (m_currentProjection)*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection - 1]->size(), m_pnProjectedPoints[m_currentProjection - 1]->normalData());PRINT_OPENGL_ERROR();
	}
	double elapsed = std::chrono::duration<double>(stop - start).count() * 1000.0;
	std::cout << "Elapsed time: " << elapsed << "ms" << std::endl;
	return elapsed;
}

double scene::projectPointsOnCPU(bool updateBuffer)
{
	//Create object for projected points
	if (m_pnProjectedPoints[m_currentProjection]->size() == 0)
	{
		delete m_pnProjectedPoints[m_currentProjection];
		if (m_random)
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject);
		else
			m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject2);
	}
	CPUTimer t;
	t.start();
	{
#pragma omp parallel for
		for (int i=0; i<(int)m_pnProjectedPoints[m_currentProjection]->validSize(); i++)
		{
			m_apss->projectCPU((*m_pnProjectedPoints[m_currentProjection])[i], (*m_pnProjectedPoints[m_currentProjection])[i], m_pnProjectedPoints[m_currentProjection]->normal(i), 1, *m_weightKernel);
		}
	}
	double elapsedTime = t.printElapsed("CPU time : ");

#ifdef CHECK
		unsigned nbOfNan;
		checkForInvalidPoints(nbOfNan);
		if (nbOfNan != 0)
		{
			cerr << "--------------------------" << endl;
			cerr << "Warning, " << nbOfNan << " NaN values found" << endl;
			cerr << "This can happen when constant shift = 0" << endl;
			cerr << "Or if the first sigmas are too small" << endl;
			cerr << "--------------------------" << endl;
		}
#endif

	if (updateBuffer)
	{
		glBindVertexArray(m_vao3);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_currentProjection * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + m_currentProjection*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->normalData());PRINT_OPENGL_ERROR();
	}
	else if (m_currentProjection != 0) //Put last pointset in buffer
	{
		glBindVertexArray(m_vao3);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_currentProjection) * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection - 1]->size(), m_pnProjectedPoints[m_currentProjection - 1]->positionData());PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + (m_currentProjection)*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection - 1]->size(), m_pnProjectedPoints[m_currentProjection - 1]->normalData());PRINT_OPENGL_ERROR();
	}
	return elapsedTime;
}

void scene::eraseFromGPU()
{
	m_apss->eraseFromGPU();
}

void scene::setNewKernel(Kernel const & newKernel)
{
	delete(m_weightKernel);
	m_weightKernel = new Kernel(newKernel);
	m_weightKernel->sigma *= m_apssOctree->getBoundingBox().radius();
	for (unsigned int i=0; i<m_weightKernel->sigmas.size(); i++)
		m_weightKernel->sigmas[i] *= m_apssOctree->getBoundingBox().radius();
	m_apss->updateKernel(*m_weightKernel);
}

double scene::resetPoints(bool useAdaptiveOctree)
{
	if (m_useAdaptiveOctree != useAdaptiveOctree)
	{
		m_useAdaptiveOctree = useAdaptiveOctree;
		return resetOctree(m_useAdaptiveOctree);
	}
	initUpsampling();
	m_apss->erasePointsFromGPU();
	if (m_useAdaptiveOctree)
	{
		Pn newPn(generateFromAdaptiveOctree());
		m_pnToProject2 = newPn;
	}
	else
	{
		Pn newPn(generateFromUsualOctree());
		m_pnToProject2 = newPn;
	}
	double mortonTime = m_pnToProject2.order();

	m_apss->erasePointsFromGPU();
	delete m_pnProjectedPoints[m_currentProjection];

	m_apss->copyPointsToGPU(m_pnToProject2.size(), m_pnToProject2.positionData());
	m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject2);

	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_currentProjection * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + m_currentProjection*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->normalData());PRINT_OPENGL_ERROR();
	return mortonTime;
}

double scene::resetOctree(bool useAdaptiveOctree)
{
	m_useAdaptiveOctree = useAdaptiveOctree;
	initUpsampling();
	if(m_useAdaptiveOctree)
		generate_adaptive_octree();
	else
		generate_usual_octree();
	return resetPoints();
}

void scene::resamplePointSet(unsigned int nbOfPoints)
{
	m_apss->eraseFromGPU();
	m_pn = m_originalPn;
	m_pn.resample(nbOfPoints);
	m_apssOctree->build(m_pn.getPositions(), m_pn.getNormals());
	m_apss->copyToGPU(*m_weightKernel, m_nbOfPoints);
}

void scene::selectPtsCloseToPlane(glm::vec3 normal, glm::vec3 pt, float threshold)
{
	for (unsigned int i=0; i<m_pn.size(); i++)
	{
		float dist = fabs(glm::dot(m_pn[i], normal) + glm::dot(pt, normal));
		if (dist < threshold)
		{
			m_selection.add(m_pn[i], m_pn.normal(i));
		}
	}
	if (m_selection.size() != 0)
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
		glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_selection.size(), m_selection.positionData());PRINT_OPENGL_ERROR();
		glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size(), 3*sizeof(float)*m_selection.size(), m_selection.normalData());PRINT_OPENGL_ERROR();
		m_pwidget->selectionNotEmpty();
	}
}

void scene::selectClosestPoint(glm::vec3 posCam, glm::vec3 farPosCam, float threshold)
{
	glm::vec3 dir = glm::normalize(farPosCam - posCam);
	glm::vec3 candidate = farPosCam;
	unsigned vector_index, pt_nber;
	for (unsigned i=0; i<m_pnProjectedPoints.size(); i++)
	{
		for (unsigned j=0; j<m_pnProjectedPoints[i]->size(); j++)
		{
			glm::vec3 p = (*m_pnProjectedPoints[i])[j];
			if (glm::dot(dir, farPosCam - p) < 0)
				continue;
			float dist = glm::length((posCam - p) - glm::dot((posCam - p), dir) * dir);
			if (dist < threshold)
			{
				if (glm::distance(p, posCam) < glm::distance(candidate, posCam))
				{
					candidate = p;
					vector_index = i;
					pt_nber = j;
				}
			}
		}
	}
	if (candidate == farPosCam)
	{
		cout << "No points were selected" << endl;
		m_selection.clear();
		m_pwidget->selectionEmpty();
		renderSelection = false;
		return;
	}
	m_selection.clear();
	m_selection.add(candidate, m_pnProjectedPoints[vector_index]->normal(pt_nber));
	m_pwidget->selectionNotEmpty();
	renderSelection = true;
	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_selection.size(), m_selection.positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size(), 3*sizeof(float)*m_selection.size(), m_selection.normalData());PRINT_OPENGL_ERROR();
}

void scene::computeNoise()
{
	m_pwidget->selectionNotEmpty();
	renderSelection = true;
	m_selection = m_pn;
	float umax = 60 * M_PI / 180;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<float> dist(0, 1);
	for (unsigned int i=0; i<m_selection.size(); i++)
	{
		const glm::vec3& e3 = m_selection.normal(i);
		glm::vec3 e1;
		if (glm::dot(e3, glm::vec3(0, 0, 1)) != 0)
			e1 = glm::normalize(glm::cross(e3, glm::vec3(0, 0, 1)));
		else
			e1 = glm::normalize(glm::cross(e3, glm::vec3(0, 1, 0)));
		const glm::vec3 e2 = glm::normalize(glm::cross(e3, e1));

		float randTheta = dist(mt) * 2 * M_PI;
		glm::vec3 randVector = glm::normalize(cos(randTheta) * e1 + sin(randTheta) * e2);
		float phiRand = umax * sqrt(dist(mt));
		m_selection.normal(i) = glm::normalize(cos(phiRand) * e3 + sin(phiRand) * randVector);
	}
	glBindBuffer(GL_ARRAY_BUFFER, m_vboSelection);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_selection.size(), m_selection.positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pn.size(), 3*sizeof(float)*m_selection.size(), m_selection.normalData());PRINT_OPENGL_ERROR();
}

void scene::changeExtension(string & path, const string & newExtension)
{
	std::size_t endPath = path.find_last_of('.');
	path = path.substr(0, endPath) + newExtension;
}

bool scene::loadParameters(string filename)
{
	changeExtension(filename, ".params");
	fstream file;
	file.open(filename, ios_base::in | ios_base::binary);
	if (!file.is_open()) //No previous parameters
		return false;

	glm::fquat q;
	file.read(reinterpret_cast<char *>(&q), sizeof(glm::fquat));
	m_pwidget->cam.setQuaternion(q);
	glm::vec3 t;
	file.read(reinterpret_cast<char *>(&t), sizeof(glm::vec3));
	m_pwidget->cam.setTranslation(t);
	float d;
	file.read(reinterpret_cast<char *>(&d), sizeof(float));
	m_pwidget->cam.setDist(d);

	float val;
	file.read(reinterpret_cast<char *>(&m_weightKernel->sigma), sizeof(float));
	file.read(reinterpret_cast<char *>(&val), sizeof(float));
	m_weightKernel->fixedSigmas.clear();
	m_weightKernel->fixedSigmas.push_back(val);
	file.read(reinterpret_cast<char *>(&val), sizeof(float));
	m_weightKernel->fixedSigmas.push_back(val);
	file.read(reinterpret_cast<char *>(&val), sizeof(float));
	m_weightKernel->fixedSigmas.push_back(val);
	file.read(reinterpret_cast<char *>(&val), sizeof(float));
	m_weightKernel->fixedSigmas.push_back(val);
	file.read(reinterpret_cast<char *>(&val), sizeof(float));
	m_weightKernel->fixedSigmas.push_back(val);
	file.read(reinterpret_cast<char *>(&m_weightKernel->beginningSigma), sizeof(float));
	file.read(reinterpret_cast<char *>(&m_weightKernel->factor), sizeof(float));
	file.read(reinterpret_cast<char *>(&m_weightKernel->constant_shift), sizeof(float));
	file.read(reinterpret_cast<char *>(&m_weightKernel->exponent), sizeof(float));
	file.read(reinterpret_cast<char *>(&val), sizeof(float));

	file.read(reinterpret_cast<char *>(&m_weightKernel->nbOfComputedSigmas), sizeof(unsigned));
	file.read(reinterpret_cast<char *>(&maxNbIterations), sizeof(unsigned));
	file.read(reinterpret_cast<char *>(&m_maxDepth), sizeof(unsigned));

	file.read(reinterpret_cast<char *>(&m_minDepth), sizeof(unsigned));
	unsigned mode;
	file.read(reinterpret_cast<char *>(&mode), sizeof(unsigned));
	m_method = (mode);
	if (m_method)
		m_apss = new Fastapss(m_apssOctree);
	else
		m_apss = new GlobalAPSS(m_apssOctree);
	m_apss->setScalingProtectionSphere(val);
	m_apss->setMinDepth(m_minDepth);

	file.read(reinterpret_cast<char *>(&m_weightKernel->type), sizeof(Kernel::KernelType));
	file.read(reinterpret_cast<char *>(&computedSigma), sizeof(bool));
	file.read(reinterpret_cast<char *>(&m_cut), sizeof(bool));

	file.read(reinterpret_cast<char *>(&m_nbOfPoints), sizeof(unsigned));
	file.read(reinterpret_cast<char *>(&m_maxUpsamples), sizeof(unsigned));
	file.read(reinterpret_cast<char *>(m_offset.data()), 6 * sizeof(float));

	file.close();

	cout << "Parameters loaded from file " << filename << endl;

	m_weightKernel->updateSigmas(computedSigma, float(m_apssOctree->getBoundingBox().radius()));
	m_pwidget->updateUI();
	return true;
}

double scene::computeDualContouring(string filename, bool recordOctree)
{
	CPUTimer t;

	cout << "Maximum depth for the upsampling octree " << m_maxDepth << endl;
	const aabb_octree octzone =
	{
		Eigen::AlignedBox3f(Eigen::Vector3f( -1.0f, -1.0f, -1.0f), Eigen::Vector3f(1.0f, 1.0f, 1.0f))
	};
	const Eigen::AlignedBox3f B = octzone.bounds();
	Eigen::Matrix4f unit2bbox = Eigen::Matrix4f::Zero();
	unit2bbox(0, 0) = B.diagonal()[0];
	unit2bbox(1, 1) = B.diagonal()[1];
	unit2bbox(2, 2) = B.diagonal()[2];
	unit2bbox(0, 3) = B.min()[0];
	unit2bbox(1, 3) = B.min()[1];
	unit2bbox(2, 3) = B.min()[2];
	unit2bbox(3, 3) = 1;

	std::vector<oriented_point> verticesAndNormals;
	std::vector<uint32_t> indices;
	cout << "Uniformizing tree ..." << endl;

	sorted_octree<char> totalUpsamplingOctree;
	glm::vec3 new_diag(1), new_origin(0);
	t.start();
	if(m_onUpsampleOctree)
	{
		OctreeGenerator octree_gen(m_apss, *m_weightKernel);
		if (m_pnProjectedPoints.size() < 1)
		{
			cerr << "Error: no projections" << endl;
			return -1;
		}
		Pn *allProjPts = new Pn(*m_pnProjectedPoints[0]);
		for (unsigned i=1; i<m_pnProjectedPoints.size(); i++)
		{
			allProjPts->add(m_pnProjectedPoints[i]->getPositions(), m_pnProjectedPoints[i]->getNormals());
		}
		std::vector<octree_id> upsampleOctree = octree_gen.generateFromPoints(m_maxDepth, *allProjPts);
		uniformizeFromPoints(upsampleOctree, totalUpsamplingOctree, new_diag, new_origin);
	}
	else
	{
		if (m_useAdaptiveOctree)
			uniformize(m_maxDepth, m_totalUpsamplingNodes, totalUpsamplingOctree, new_diag, new_origin, m_cut);
		else
			uniformizeFromPoints(m_usualOctree, totalUpsamplingOctree, new_diag, new_origin);
	}
	m_uniformizeTime = t.printElapsedAndReset("[ OK ] ");
	cout << totalUpsamplingOctree.size() << " nodes" << endl;

	Eigen::Matrix4f sortedOctreeToUnit = Eigen::Matrix4f::Zero();
	sortedOctreeToUnit(0, 0) = new_diag[0];
	sortedOctreeToUnit(1, 1) = new_diag[1];
	sortedOctreeToUnit(2, 2) = new_diag[2];
	sortedOctreeToUnit(0, 3) = new_origin[0];
	sortedOctreeToUnit(1, 3) = new_origin[1];
	sortedOctreeToUnit(2, 3) = new_origin[2];
	sortedOctreeToUnit(3, 3) = 1;

	Eigen::Matrix4f new2bbox = unit2bbox * sortedOctreeToUnit;

	if (recordOctree)
	{
		cout << "Saving octree ..." << endl;
		Eigen::Matrix4f bbox2orig = Eigen::Matrix4f::Zero();
		bbox2orig(0, 0) = m_originalPn.getNorme();
		bbox2orig(1, 1) = m_originalPn.getNorme();
		bbox2orig(2, 2) = m_originalPn.getNorme();
		bbox2orig(0, 3) = m_originalPn.getOffset().x;
		bbox2orig(1, 3) = m_originalPn.getOffset().y;
		bbox2orig(2, 3) = m_originalPn.getOffset().z;
		bbox2orig(3, 3) = 1;
		const Eigen::Matrix4f transform = bbox2orig * new2bbox;
		std::vector<octree_id> N;
		for (unsigned j=0; j <= unsigned(totalUpsamplingOctree.max_depth()); j++)
		{
			N.clear();
			std::cout << "depth: " << j << "...";
			for(const octree_id& n : totalUpsamplingOctree.all_nodes(j))
				N.push_back(n);
			std::cout << " (" << N.size() << " nodes)... ";
			string filenameO = filename;
			changeExtension(filenameO, "Octree" + to_string(j) + ".ply");
			save_octree_ply(filenameO, N, node_color_by_depth, transform);
			std::cout << " saved\n";
		}
		string filenameL = filename;
		changeExtension(filenameL, "Leaves.ply");
		save_octree_ply(filenameL, m_adaptiveOctree, node_color_by_depth, bbox2orig * unit2bbox);
		cout << "Files " << "OctreeX.ply" << " and " << filenameL << " have been written" << endl;
	}

	auto surface = [this, &new2bbox](const Eigen::Vector3f& p) -> float
	{
		glm::vec3 pproj, nproj;
		Eigen::Vector4f pt = new2bbox * Eigen::Vector4f(p.x(), p.y(), p.z(), 1.f);//We go to the frame of the Pn
		glm::vec3 originPt = glm::vec3(pt.x(), pt.y(), pt.z());
		m_apss->projectCPU(originPt, pproj, nproj, 1, *m_weightKernel);
		return glm::dot(( originPt - pproj ), nproj);
	};

	auto computeNormals = [](const Eigen::Vector3f& p) -> Eigen::Vector3f
	{
		return Eigen::Vector3f(0, 0, 0);
	};

	//Precompute corners projections in parallel using CUDA for better performances with the dual contouring
	std::unordered_map<Eigen::Vector3i, float, vec3i_hash> corners;
	size_t octSize = totalUpsamplingOctree.size();
	std::vector<glm::vec3> listOfCorners(octSize * 8), listOfNormals(octSize * 8), listOfProjectedCorners(octSize * 8);
	std::vector<Eigen::Vector3i> allCorners(octSize * 8);
#pragma omp parallel for
	for (unsigned i=0; i<octSize; i++)
	{
		const octree_id n = totalUpsamplingOctree.id(i);
		for (unsigned c=0; c<8; c++)
		{
			const float w = 1.0f / (1 << totalUpsamplingOctree.max_depth());
			const int delta = 1 << (totalUpsamplingOctree.max_depth() - n.depth);
			const Eigen::Vector3i p = delta * n.corner(c);
			allCorners[i * 8 + c] = p;
			const Eigen::Vector3f pf = p.cast<float>() * w;
			const Eigen::Vector4f pt = new2bbox * Eigen::Vector4f(pf.x(), pf.y(), pf.z(), 1.f);
			listOfCorners[i * 8 + c] = glm::vec3(pt.x(), pt.y(), pt.z());
		}
	}
	cout << "Number of corners to project: " << listOfCorners.size() << endl;

	{
		m_apss->erasePointsFromGPU();
		m_apss->copyPointsToGPU(listOfCorners.size(), listOfCorners.data());
		m_apss->project(listOfCorners.size(), listOfProjectedCorners.data(), listOfNormals.data(), 1, *m_weightKernel);
	}

	for (unsigned i=0; i<listOfCorners.size(); i++)
	{
		if (glm::length2(listOfNormals[i]) < 0.5f)//Normal is (0, 0, 0) if the projection failed
			corners[allCorners[i]] = std::nanf("");
		else
			corners[allCorners[i]] = glm::dot(( glm::vec3(listOfCorners[i]) - listOfProjectedCorners[i] ), listOfNormals[i]);
	}

	cout << "Computing dual contouring ..." << endl;
	t.start();
	dual_contouring(totalUpsamplingOctree, surface, computeNormals, corners).extract(verticesAndNormals, indices);
	t.stop();
	cout << "     Done" << endl;

	std::vector<glm::vec3> vertices(verticesAndNormals.size());
	std::vector<glm::vec3> normals(verticesAndNormals.size());
#pragma omp parallel for
	for (unsigned i=0; i<verticesAndNormals.size(); i++)
	{
		Eigen::Vector4f v = new2bbox * Eigen::Vector4f(verticesAndNormals[i].position.x(), verticesAndNormals[i].position.y(), verticesAndNormals[i].position.z(), 1.f);
		vertices[i] = glm::vec3(v.x(), v.y(), v.z());
		normals[i] = glm::vec3(verticesAndNormals[i].normal.x(), verticesAndNormals[i].normal.y(), verticesAndNormals[i].normal.z());
		if (isnan(v.x()) || isnan(v.y()) || isnan(v.z()))
			cerr << "NAN" << endl;
	}

	cout << vertices.size() << " points to project ..." << endl;

	//Projection
	{
		if (m_apssOctree->isAPSS())
		{
#pragma omp parallel for
			for (unsigned i=0; i<vertices.size(); i++)
				m_apss->projectCPU(vertices[i], vertices[i], normals[i], 3, *m_weightKernel);
		}
		else
		{
			m_apss->erasePointsFromGPU();
			m_apss->copyPointsToGPU(vertices.size(), vertices.data());
			m_apss->project(vertices.size(), vertices.data(), normals.data(), 3, *m_weightKernel);
		}
	}
	//Smoothing and Subdivision
	for (unsigned i=0; i<2; i++)
	{
		mesh::smoothMesh(vertices, indices);
		mesh::subdivideMesh(vertices, indices);
		normals.resize(vertices.size());
		//Projection
		{
			if (m_apssOctree->isAPSS())
			{
#pragma omp parallel for
				for (unsigned i=0; i<vertices.size(); i++)
					m_apss->projectCPU(vertices[i], vertices[i], normals[i], 3, *m_weightKernel);
			}
			else
			{
				m_apss->erasePointsFromGPU();
				m_apss->copyPointsToGPU(vertices.size(), vertices.data());
				m_apss->project(vertices.size(), vertices.data(), normals.data(), 3, *m_weightKernel);
			}
		}
	}
	if (m_currentProjection + 1 < m_maxUpsamples && m_pnProjectedPoints[m_currentProjection]->size() != 0)
	{
		m_apss->erasePointsFromGPU();
		m_apss->copyPointsToGPU(m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());

	}
	else if (m_currentProjection == 0)
	{
		m_apss->erasePointsFromGPU();
		m_apss->copyPointsToGPU(m_pnToProject2.size(), m_pnToProject2.positionData());
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_vboDualContouring);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*vertices.size()*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	glBindVertexArray(m_vaoDualContouring);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboDualContouring);
	GLint loc12=glGetAttribLocation(m_triangleProgram, "position");
	glVertexAttribPointer(loc12, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc12);
	GLint loc13=glGetAttribLocation(m_triangleProgram, "normal");
	glVertexAttribPointer(loc13, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*vertices.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc13);

	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*vertices.size(), vertices.data()); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*vertices.size(), 3*sizeof(float)*vertices.size(), normals.data()); PRINT_OPENGL_ERROR();

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboDualContouring);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(uint32_t), indices.data(), GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	m_nbOfIndices = indices.size();
	renderDualContour = true;

	fstream file;
	if (filename == "")
		filename = m_pwidget->getFilename();
	string filenameDC = filename;
	file.open(filenameDC, ios_base::out | ios_base::binary);
	if (!file.is_open())
	{
		std::cerr << "Impossible to write file " << filenameDC << std::endl;
		return t.printElapsed("Time for dual contouring");
	}
	file << "ply" << endl;
	file << "format binary_little_endian 1.0" << endl;
	file << "comment Created by FastAPSS" << endl;
	file << "comment Kernel parameters:" << endl;
	switch (m_weightKernel->type) {
	case Kernel::GAUSSIAN_MULTIPLE:
		file << "comment Multiple Gaussian" << endl;
		file << "comment First sigma  " << m_weightKernel->beginningSigma << endl;
		file << "comment Factor " << m_weightKernel->factor << endl;
		file << "comment Number of sigmas " << m_weightKernel->sigmas.size() << endl;
		file << "comment Constant_shift (unused) " << m_weightKernel->constant_shift << endl;
		break;
	case Kernel::SINGULAR:
		file << "comment Singular" << endl;
		file << "comment Exponent " << m_weightKernel->exponent << endl;
		file << "comment Constant_shift " << m_weightKernel->constant_shift << endl;
		break;
	case Kernel::GAUSSIAN:
		file << "comment Gaussian" << endl;
		file << "comment Sigma " << m_weightKernel->sigma / m_apssOctree->getBoundingBox().radius() << endl;
		file << "comment Constant_shift " << m_weightKernel->constant_shift << endl;
		break;
	default:
		file << "comment Kernel unknown" << endl;
		break;
	}
	file << "comment Octree:" << endl;
	if (m_onUpsampleOctree)
		file << "comment Built on upsampled points" << endl;
	else
	{
		if (m_useAdaptiveOctree)
			file << "comment Built using the adaptive octree method" << endl;
		else
			file << "comment Built on input points" << endl;
	}
	file << "element vertex " << int(vertices.size()) << endl;
	file << "property float x" << endl;
	file << "property float y" << endl;
	file << "property float z" << endl;
	file << "property float nx" << endl;
	file << "property float ny" << endl;
	file << "property float nz" << endl;
	file << "element face " << int(indices.size() / 3) << endl;
	file << "property list uchar int vertex_indices" << endl;
	file << "end_header" << endl;
	for (unsigned int i=0; i<vertices.size(); i++)
	{
		glm::vec3 realVertex = (vertices[i] * m_originalPn.getNorme()) + m_originalPn.getOffset();
		file.write(reinterpret_cast<const char *>(&(realVertex)), sizeof (glm::vec3));
		file.write(reinterpret_cast<const char *>(&(normals[i])), sizeof (glm::vec3));
	}
	int triangleSize = 3;
	for (unsigned int i=0; i<indices.size()/3; i++)
	{
		file.write(reinterpret_cast<const char *>(&triangleSize), sizeof (unsigned char));
		file.write(reinterpret_cast<const char *>(&(indices[i*3])), sizeof (int));
		file.write(reinterpret_cast<const char *>(&(indices[i*3 + 1])), sizeof (int));
		file.write(reinterpret_cast<const char *>(&(indices[i*3 + 2])), sizeof (int));
	}
	file.close();
	cout << "File " << filenameDC << " saved" << endl;
	return t.printElapsed("Time for dual contouring");
}

void scene::createTrianglesForMPU()
{
	vector<glm::vec3> vertices(m_originalPn.size() * 3);

	float f = 0.001f;

#pragma omp parallel for
	for (int i=0; i<(int)m_originalPn.size(); i++)
	{
		glm::vec3 a(0, 0, 1);
		if (abs(glm::dot(a, m_originalPn.normal(i))) < 0.00001f)
			a = glm::vec3(0, 1, 0);
		vertices[i * 3] = m_originalPn[i] + glm::normalize(glm::cross(a, m_originalPn.normal(i))) * f;
		glm::mat4 r = glm::rotate(float(2.f * M_PI / 3.f), m_originalPn.normal(i));
		vertices[i * 3 + 1] = m_originalPn[i] + glm::vec3(glm::vec4(vertices[i * 3] - m_originalPn[i], 1) * r);
		vertices[i * 3 + 2] = m_originalPn[i] + glm::vec3(glm::vec4(vertices[i * 3 + 1] - m_originalPn[i], 1) * r);
	}

	fstream file;
	string filename = m_pwidget->getFilename();
	changeExtension(filename, "Triangles.ply");
	file.open(filename, ios_base::out | ios_base::binary);
	if (!file.is_open())
	{
		std::cerr << "Impossible to write file " << filename << std::endl;
		return;
	}
	file << "ply" << endl;
	file << "format binary_little_endian 1.0" << endl;
	file << "element vertex " << int(vertices.size()) << endl;
	file << "property float x" << endl;
	file << "property float y" << endl;
	file << "property float z" << endl;
	file << "property float nx" << endl;
	file << "property float ny" << endl;
	file << "property float nz" << endl;
	file << "element face " << int(vertices.size() / 3) << endl;
	file << "property list uchar int vertex_indices" << endl;
	file << "end_header" << endl;
	for (unsigned int i=0; i<vertices.size(); i++)
	{
		glm::vec3 realVertex = (vertices[i] * m_originalPn.getNorme()) + m_originalPn.getOffset();
		file.write(reinterpret_cast<const char *>(&(realVertex)), sizeof (glm::vec3));
		file.write(reinterpret_cast<const char *>(&(m_originalPn.normal(i / 3))), sizeof (glm::vec3));
	}
	int triangleSize = 3;
	int v;
	for (int i=0; i<int(vertices.size())/3; i++)
	{
		file.write(reinterpret_cast<const char *>(&triangleSize), sizeof (unsigned char));
		v = i*3;
		file.write(reinterpret_cast<const char *>(&v), sizeof (int));
		v = i*3+1;
		file.write(reinterpret_cast<const char *>(&v), sizeof (int));
		v = i*3+2;
		file.write(reinterpret_cast<const char *>(&v), sizeof (int));
	}
	file.close();
	cout << "File " << filename << " saved" << endl;
}

void scene::saveProjectedPoints(string filename)
{
	if (m_pnProjectedPoints.size() == 0)
	{
		cerr << "No projected points to record" << endl;
		return;
	}
	unsigned nbOfBatches = std::min(m_currentProjection + 1, m_maxUpsamples);
	//Only ply files here
	fstream file;
	file.open(filename, ios_base::out | ios_base::binary);
	if (!file.is_open())
	{
		std::cerr << "Impossible to write file " << filename << std::endl;
		return;
	}
	file << "ply" << endl;
	file << "format binary_little_endian 1.0" << endl;
	int nbOfPoints = int(nbOfBatches * m_pnProjectedPoints[0]->size());
	file << "element vertex " << nbOfPoints << endl;
	file << "property float x" << endl;
	file << "property float y" << endl;
	file << "property float z" << endl;
	file << "property float nx" << endl;
	file << "property float ny" << endl;
	file << "property float nz" << endl;
	file << "end_header" << endl;
	for (unsigned int i=0; i<nbOfBatches; i++)
	{
		for (unsigned int j=0; j<m_pnProjectedPoints[i]->size(); j++)
		{
			glm::vec3 position = ((*m_pnProjectedPoints[i])[j] * m_originalPn.getNorme()) + m_originalPn.getOffset();
			file.write(reinterpret_cast<const char *>(&(position)), sizeof (glm::vec3));
			file.write(reinterpret_cast<const char *>(&(m_pnProjectedPoints[i]->normal(j))), sizeof (glm::vec3));
		}
	}
	file.close();

	cout << "File " << filename << " has been saved !" << endl;
}

void scene::comparePointSet(string filename)
{
	m_pnForComparison.load(filename, m_pwidget->getRadiusFactor(), false, m_pn.getOffset(), m_pn.getNorme());


	Pn projectedPoints = m_pnForComparison;
	m_apss->project(projectedPoints.size(), projectedPoints.positionData(), projectedPoints.normalData(), 1, *m_weightKernel);

	vector<float> perPointValue(projectedPoints.size());
#pragma omp parallel for
	for (int i=0; i<(int)projectedPoints.size(); i++)
	{
		perPointValue[i] = glm::dot(m_pnForComparison[i] - projectedPoints[i], projectedPoints.normal(i));
	}

	float min = perPointValue[0];
	float max = perPointValue[0];
	for (unsigned i=0; i<projectedPoints.size(); i++)
	{
		if (perPointValue[i] < min)
			min = perPointValue[i];
		if (perPointValue[i] > max)
			max = perPointValue[i];
	}
	float norme = std::max(abs(max), abs(min));
#pragma omp parallel for
	for (int i=0; i<(int)projectedPoints.size(); i++)
	{
		perPointValue[i] = (perPointValue[i]) / norme;
	}

	glBindVertexArray(m_vaoComparison);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboComparison);
	glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnForComparison.size()*3, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	GLint loc1=glGetAttribLocation(m_perPointColorProgram, "position");
	glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc1);
	GLint loc2=glGetAttribLocation(m_perPointColorProgram, "normal");
	glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_pnForComparison.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc2);
	GLint loc3=glGetAttribLocation(m_perPointColorProgram, "perPointColor");
	glVertexAttribPointer(loc3, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(6*sizeof(float)*m_pnForComparison.size()));PRINT_OPENGL_ERROR();
	glEnableVertexAttribArray(loc3);

	glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*m_pnForComparison.size(), m_pnForComparison.positionData()); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_pnForComparison.size(), 3*sizeof(float)*m_pnForComparison.size(), m_pnForComparison.normalData()); PRINT_OPENGL_ERROR();

	vector<glm::vec3> colors(m_pnForComparison.size());
#pragma omp parallel for
	for (int i=0; i<(int)m_pnForComparison.size(); i++)
	{
		colors[i] = Pn::GetColour(perPointValue[i], -1, 1);
	}
	glBufferSubData(GL_ARRAY_BUFFER, 6*sizeof(float)*m_pnForComparison.size(), 3*sizeof(float)*m_pnForComparison.size(), colors.data()); PRINT_OPENGL_ERROR();

	renderComparison = true;
}

void scene::checkForInvalidPoints(unsigned & nbOfInvalidPts)
{
	nbOfInvalidPts = 0;
	unsigned i = m_currentProjection;
	std::vector<unsigned> invalidPts;
	for (unsigned j=0; j<m_pnProjectedPoints[i]->validSize(); j++)
	{
		bool is_invalid = true;
		for (int k=0; k<3; k++)
			is_invalid = is_invalid && (m_pnProjectedPoints[i]->normal(j)[k] == 0);
		if (is_invalid)
		{
			nbOfInvalidPts++;
			invalidPts.push_back(j);
		}
	}
	//Invalid points are put at the end of the vector, so that they are not used any more
#pragma omp parallel for
	for (unsigned j=0; j<invalidPts.size(); j++)
	{
		const glm::vec3 temp = (*m_pnProjectedPoints[i])[invalidPts[j]];
		(*m_pnProjectedPoints[i])[invalidPts[j]] = (*m_pnProjectedPoints[i])[m_pnProjectedPoints[i]->validSize() - 1 - j];
		(*m_pnProjectedPoints[i])[m_pnProjectedPoints[i]->validSize() - 1 - j] = temp;
		const glm::vec3 temp2 = m_pnProjectedPoints[i]->normal(invalidPts[j]);
		m_pnProjectedPoints[i]->normal(invalidPts[j]) = m_pnProjectedPoints[i]->normal(m_pnProjectedPoints[i]->validSize() - 1 - j);
		m_pnProjectedPoints[i]->normal(m_pnProjectedPoints[i]->validSize() - 1 - j) = temp2;
	}
	if (invalidPts.size() > 0)
		m_apss->reorganizePoints(&invalidPts, m_pnProjectedPoints[i]->validSize());
	m_pnProjectedPoints[i]->setValidSize(m_pnProjectedPoints[i]->validSize() - nbOfInvalidPts);
	m_nbOfInvalidPts += nbOfInvalidPts;
}

void scene::loadPointsToProject()
{
	m_pnToProject.load("PointsToProject.ply", m_pwidget->getRadiusFactor(), false);
	initUpsampling();
	delete m_pnProjectedPoints[m_currentProjection];
	m_pnProjectedPoints[m_currentProjection] = new Pn(m_pnToProject);
	m_apss->erasePointsFromGPU();
	m_apss->copyPointsToGPU(m_pnToProject.size(), m_pnToProject.positionData());
	glBindBuffer(GL_ARRAY_BUFFER, m_vboPointSetProjected);                 PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_currentProjection * m_pnToProject.size(), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->positionData());PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*(m_pnProjectedPoints[m_currentProjection]->size()*m_maxUpsamples + m_currentProjection*m_pnToProject.size()), 3*sizeof(float)*m_pnProjectedPoints[m_currentProjection]->size(), m_pnProjectedPoints[m_currentProjection]->normalData());PRINT_OPENGL_ERROR();
}

void scene::mortonTest()
{
	cout << "Beginning of Morton test..." << endl;
	unsigned maxUpsamples = m_maxUpsamples;
	unsigned nbOfPts = m_nbOfPoints;
	m_maxUpsamples = 1;
	fstream file;
	file.open("Morton.txt", ios_base::out);
	for (unsigned i=100; i<nbOfPts+1; i+=100)
	{
		cout << "Current number of points: " << i << endl;
		m_nbOfPoints = i;
		file << i << " " << updateBuffers(true) << " "; //First, the number of points, then the time for ordering with the Morton Code
		double meanTime = 0.0;
		for (unsigned j=0; j<10; j++)
			meanTime += projectPoints(true);
		meanTime /= 10.0;
		file << meanTime << " "; //Projection time when points are ordered

		updateBuffers(false);
		meanTime = 0.0;
		for (unsigned j=0; j<10; j++)
			meanTime += projectPoints(true);
		meanTime /= 10.0;
		file << meanTime << endl; //Projection time when points are not ordered
	}
	file.close();
	m_maxUpsamples = maxUpsamples;
	m_nbOfPoints = nbOfPts;
	cout << "Morton test finished" << endl;
}
