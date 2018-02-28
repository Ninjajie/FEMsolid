# ifndef FEMSolid_H
# define FEMSolid_H

# define TETLIBRARY

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
class FEMSolidSolver
{
private:
	fReal timeStep;
	fReal framePeriod;

	long long steps;
	long long stepsPerFrame;

	std::vector<Eigen::Vector4i> tetraIndices;

	std::vector<vec3> positions;
	std::vector<vec3> velocities;
	std::vector<vec3> elasticForces;
	std::vector<vec3> bodyForces;
	vec3 sphereOrigin;
	vec3 sphereVelocity;
	fReal sphereMass;
	std::vector<fReal> masses;

	std::vector<mat3> Dm;
	std::vector<mat3> Bm;
	std::vector<fReal> We;

	void preCompute();

	void computeElasticForce();

	void computeBodyForce();

	void save2File(std::string path);

	void sphereBound(fReal idx, vec3 origin);

	void save2filesSphere(vec3 o, fReal r, std::string path);
public:
	//以后把这个弄个虚函数就好了
	virtual mat3 computeP(mat3 dst);

	void stepForward();

	static FEMSolidSolver* createFromObj(const std::string& filename, fReal timeStep, fReal framePeriod);

	static FEMSolidSolver* createFromCube(fReal timeStep, fReal framePeriod);

	static FEMSolidSolver* createForDebugging(fReal timeStep, fReal framePeriod);

	FEMSolidSolver(fReal timeStep, fReal framePeriod);

	long long getCurrentIterations();
};

# endif