# include "FEMSolid.h"

const fReal size = 0.3;

FEMSolidSolver* FEMSolidSolver::createForDebugging(fReal timeStep, fReal framePeriod)
{
	tetgenio p;
	FEMSolidSolver* solver = new FEMSolidSolver(p, timeStep, framePeriod);
	solver->preAllocate(5, 8);

	solver->tetraIndices[0] = Eigen::Vector4i(0, 1, 3, 4);
	solver->tetraIndices[1] = Eigen::Vector4i(1, 4, 5, 6);
	solver->tetraIndices[2] = Eigen::Vector4i(1, 2, 3, 6);
	solver->tetraIndices[3] = Eigen::Vector4i(3, 4, 6, 7);
	solver->tetraIndices[4] = Eigen::Vector4i(1, 3, 4, 6);
	
	solver->positions[0] = vec3(size, -size, -size);
	solver->positions[1] = vec3(size, -size, size);
	solver->positions[2] = vec3(size, size, size);
	solver->positions[3] = vec3(size, size, -size);
	solver->positions[4] = vec3(-size, -size, -size);
	solver->positions[5] = vec3(-size, -size, +size);
	solver->positions[6] = vec3(-size, size, size);
	solver->positions[7] = vec3(-size, size, -size);

	for (int i = 0; i < 8; i++)
	{
		solver->velocities[i] = (zeroVec);
	}

	solver->preCompute();

	return solver;
}