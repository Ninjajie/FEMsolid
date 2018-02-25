# ifndef FEMSolid_H
# define FEMSolid_H

# define TETLIBRARY

/*Put every last 3rd party we need in this header*/
# include <iostream>
# include <cublas.h>
# include <cuda.h>
# include <Eigen\Eigen>
# include <tetgen.h>
# include <array>
# include <vector>
# include <poly.hpp>

class FEMSolidSolver
{
private:
	double timeStep;
	double framePeriod;

	long long steps;
	long long stepsPerFrame;

	std::vector<Eigen::Vector4i> tetraIndices;

	std::vector<Eigen::Vector3d> positions;
	std::vector<Eigen::Vector3d> velocities;
	std::vector<Eigen::Vector3d> elasticForces;
	std::vector<Eigen::Vector3d> bodyForces;
	std::vector<double> masses;

	std::vector<Eigen::Matrix3d> Dm;
	std::vector<Eigen::Matrix3d> Bm;
	std::vector<double> We;

	void preCompute();

	void computeElasticForce();

	void computeBodyForce();

public:
	//以后把这个弄个虚函数就好了
	virtual Eigen::Matrix3d computeP(Eigen::Matrix3d dst);

	void stepForward();

	static FEMSolidSolver* createFromObj(const std::string& filename, double timeStep, double framePeriod);

	static FEMSolidSolver* createFromCube(double timeStep, double framePeriod);

	FEMSolidSolver(double timeStep, double framePeriod);

	long long getCurrentIterations();
};

# endif