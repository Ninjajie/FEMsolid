# include "FEMSolid.h"

void tetrahedralizeCube(tetgenio& out);

FEMSolidSolver::FEMSolidSolver(double timeStep, double framePeriod)
	: timeStep(timeStep), framePeriod(framePeriod), steps(0), stepsPerFrame(static_cast<long long>(framePeriod / timeStep))
{}

FEMSolidSolver* FEMSolidSolver::createFromCube(double timeStep, double framePeriod)
{
	tetgenio cube;
	tetrahedralizeCube(cube);

	FEMSolidSolver* solver = new FEMSolidSolver(timeStep, framePeriod);

	for (int i = 0; i < cube.numberoftetrahedra; i++)
	{
		solver->tetraIndices.push_back(Eigen::Vector4i(cube.tetrahedronlist[4 * i],
													   cube.tetrahedronlist[4 * i + 1],
													   cube.tetrahedronlist[4 * i + 2],
													   cube.tetrahedronlist[4 * i + 3]));
	}

	for (int i = 0; i < cube.numberofpoints; i++)
	{
		Eigen::Vector3d zero = Eigen::Vector3d();
		solver->positions.push_back(Eigen::Vector3d(cube.pointlist[3 * i],
													cube.pointlist[3 * i + 1],
													cube.pointlist[3 * i + 2]));
		solver->velocities.push_back(zero);

		//TODO: 质量怎么办？
		solver->masses.push_back(1.0);
		//解的时候再把力填进去更好
		//solver->forces.push_back(zero);
	}

	solver->preCompute();

	return solver;
}

void FEMSolidSolver::preCompute()
{
	for (int i = 0; i < tetraIndices.size(); i++)
	{
		Eigen::Matrix3d dmt;
		Eigen::Vector3d& lastCorner = positions[tetraIndices[i][3]];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Eigen::Vector3d& corner = positions[tetraIndices[i][cornerInd]];

			dmt.row(cornerInd) << corner - lastCorner;
		}
		Dm.push_back(dmt);
		Bm.push_back(dmt.inverse());
		We.push_back(dmt.determinant());
	}
}

Eigen::Matrix3d FEMSolidSolver::computeP(Eigen::Matrix3d dst)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(dst, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();
	Eigen::Vector3d Sigma = svd.singularValues();

	const double E = 0.1;
	const double NU = 0.4999;

	double mu = E / (2 * (1 + NU));
	double lambda = (E * NU) / ((1 + NU) * (1 - 2 * NU));
	Eigen::Matrix3d R = U*V;
	Eigen::Matrix3d P = 2 * mu * (dst - R) + lambda * (R.transpose() * dst - Eigen::MatrixXd::Identity(3, 3)).trace() * R;
	return P;
}

void FEMSolidSolver::computeBodyForce()
{
	
}

void FEMSolidSolver::stepForward()
{
	Eigen::Vector3d zeroVec = Eigen::Vector3d(0.0, 0.0, 0.0);
	bodyForces.clear();
	elasticForces.clear();
	for (size_t pointInd = 0; pointInd != positions.size(); ++pointInd)
	{
		bodyForces.push_back(zeroVec);
		elasticForces.push_back(zeroVec);
	}

	this->computeBodyForce();
	this->computeElasticForce();

	for (size_t pointInd = 0; pointInd != positions.size(); ++pointInd)
	{
		Eigen::Vector3d totalForce = bodyForces[pointInd] + elasticForces[pointInd];
		velocities[pointInd] += this->timeStep * totalForce / masses[pointInd];
		positions[pointInd] += this->timeStep * velocities[pointInd];
	}

	if (this->steps % this->stepsPerFrame == 0)
	{
		//Save file
	}
	this->steps++;
}

void FEMSolidSolver::computeElasticForce()
{
	for (int tetInd = 0; tetInd < tetraIndices.size(); tetInd++)
	{
		Eigen::Matrix3d dst;
		Eigen::Vector3d& lastCorner = positions[tetraIndices[tetInd][3]];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Eigen::Vector3d& corner = positions[tetraIndices[tetInd][cornerInd]];
			dst.row(cornerInd) << corner - lastCorner;
		}
		dst *= Bm[tetInd];

		// todo: compute P
		Eigen::Matrix3d P = this->computeP(dst);

		// todo: what if rotation matrices are reflections
		if (dst.determinant() < 0)
		{

		}
		
		// populate H
		Eigen::Matrix3d H = -We[tetInd] * P * Bm[tetInd].transpose();

		Eigen::Vector3d f3 = Eigen::Vector3d();
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			elasticForces[tetraIndices[tetInd][cornerInd]] = H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		// 等等……啥？应该是+=吧
		elasticForces[tetraIndices[tetInd][3]] += f3;
	}
}