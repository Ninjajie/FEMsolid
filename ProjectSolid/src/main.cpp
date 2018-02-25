# include <iostream>
# include "FEMSolid.h"

using namespace Eigen;

int main(int argc, char** argv)
{
	FEMSolidSolver* solver = FEMSolidSolver::createFromCube(0.01f, 0.1f);

	for (long long i = 0; i != 100; ++i)
	{
		solver->stepForward();
	}

	delete solver;
}