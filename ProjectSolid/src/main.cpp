# include <iostream>
# include "FEMSolid.h"

using namespace Eigen;

int main(int argc, char** argv)
{
	FEMSolidSolver* solver = FEMSolidSolver::createForDebugging(0.0001f, 0.01f);

	for (long long i = 0; i != 10000; ++i)
	{
		solver->stepForward();
	}

	delete solver;
}