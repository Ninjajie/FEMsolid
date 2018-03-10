# ifndef FEMSolid_H
# define FEMSolid_H

# define TETLIBRARY
# define TINYOBJLOADER_IMPLEMENTATION
/// 如果不想要并行化代码就把这行注释掉
# define OMParallelize

# ifdef OMParallelize
# define TOTALThreads 16
# endif

# define GravityAcc 9.8

/*Put every last 3rd party we need in this header*/
# include <iostream>
# include <cublas.h>
# include <cuda.h>
# include <Eigen\Eigen>
# include <tetgen.h>
# include <array>
# include <vector>
# include <map>
# include <poly.hpp>
# include <omp.h>

# define DOUBLE_Float

# ifdef SINGLE_Float
typedef Eigen::Vector3f vec3;
typedef Eigen::Matrix3f mat3;
typedef float fReal;
# endif

# ifdef DOUBLE_Float
typedef Eigen::Vector3d vec3;
typedef Eigen::Matrix3d mat3;
typedef double fReal;
# endif


const fReal E = 5000;
const fReal NU = 0.3;
const fReal Density = 6.8;
const fReal PI = 3.1415926;

const vec3 zeroVec = vec3(0.0, 0.0, 0.0);


class FEMSolidSolver
{
private:
	fReal timeStep;
	fReal framePeriod;
	long long steps;
	long long stepsPerFrame;
	long long frames;


	Eigen::Vector4i* tetraIndices;
	mat3* Dm;
	mat3* Bm;
	fReal* We;
	size_t numOfTets;


	vec3* positions;
	vec3* velocities;
	vec3* elasticForces;
	vec3* bodyForces;
	fReal* masses;
	vec3 sphereOrigin;
	vec3 sphereVelocity;
	fReal sphereMass;
	size_t numOfVerts;


	size_t numOfObjTriangles;
	int* objTriangleIndices;
	std::map<int, int> objVerticesTable;


	void preAllocate(int numOfTets, int numOfVerts);

	void preCompute();

	void preFill(tetgenio& mesh);

	void setInitial();

	void computeElasticForce();

	void computeBodyForce();

	void solveForBoundary(int pointInd);

	void save2Poly(std::string path);

	void save2Obj(std::string path);

	void sphereBound(int idx, vec3 origin);

	void save2filesSphere(vec3 o, fReal r, std::string path);

	FEMSolidSolver(tetgenio& mesh, fReal timeStep, fReal framePeriod);


public:
	//以后把这个弄个虚函数就好了
	virtual mat3 computeP(mat3 dst);

	void stepForward();

	static FEMSolidSolver* createFromObj(const std::string& filename, fReal timeStep, fReal framePeriod);

	static FEMSolidSolver* createFromCube(fReal timeStep, fReal framePeriod);

	static FEMSolidSolver* createForDebugging(fReal timeStep, fReal framePeriod);

	long long getCurrentIterations();

	long long getCurrentFrames();

	~FEMSolidSolver();
};

# endif