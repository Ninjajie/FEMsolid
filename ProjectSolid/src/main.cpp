# include <iostream>
# include "FEMSolid.h"

using namespace Eigen;

int main(int argc, char** argv)
{
	//FEMSolidSolver* solver = FEMSolidSolver::createFromObj("cone.obj", 0.0001f, 0.01f);
	FEMSolidSolver* solver = FEMSolidSolver::createFromCube(0.0001f, 0.01f);

	for (long long i = 0; i != 60000; ++i)
	{
		solver->stepForward();
	}

	delete solver;
}