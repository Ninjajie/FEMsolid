# include "FEMSolid.h"

const double size = 1.0;

FEMSolidSolver* FEMSolidSolver::createForDebugging(double timeStep, double framePeriod)
{
	FEMSolidSolver* solver = new FEMSolidSolver(timeStep, framePeriod);

	solver->tetraIndices.push_back(Eigen::Vector4i(0, 1, 3, 4));
	solver->tetraIndices.push_back(Eigen::Vector4i(1, 4, 5, 6));
	solver->tetraIndices.push_back(Eigen::Vector4i(1, 2, 3, 6));
	solver->tetraIndices.push_back(Eigen::Vector4i(3, 4, 6, 7));
	solver->tetraIndices.push_back(Eigen::Vector4i(1, 3, 4, 6));
	//solver->tetraIndices.push_back(Eigen::Vector4i(1, 4, 5, 7));
	
	solver->positions.push_back(Eigen::Vector3d(size, -size, -size));
	solver->positions.push_back(Eigen::Vector3d(size, -size, size));
	solver->positions.push_back(Eigen::Vector3d(size, size, size));
	solver->positions.push_back(Eigen::Vector3d(size, size, -size));
	solver->positions.push_back(Eigen::Vector3d(-size, -size, -size));
	solver->positions.push_back(Eigen::Vector3d(-size, -size, +size));
	solver->positions.push_back(Eigen::Vector3d(-size, size, size));
	solver->positions.push_back(Eigen::Vector3d(-size, size, -size));

	for (int i = 0; i < 8; i++)
	{
		Eigen::Vector3d zero = Eigen::Vector3d();
		solver->velocities.push_back(zero);
		//TODO: 质量怎么办？
		solver->masses.push_back(1.0);
	}

	solver->preCompute();

	return solver;
}