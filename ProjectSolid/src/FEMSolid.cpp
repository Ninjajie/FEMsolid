# include "FEMSolid.h"

void tetrahedralizeCube(tetgenio& out);

FEMSolidSolver::FEMSolidSolver(double timeStep, double framePeriod)
	: timeStep(timeStep), framePeriod(framePeriod), steps(0), stepsPerFrame(static_cast<long long>(framePeriod / timeStep))
{
# ifdef OMParallelize
	omp_set_num_threads(TOTALThreads);
# endif
}

long long FEMSolidSolver::getCurrentIterations()
{
	return this->steps;
}

FEMSolidSolver* FEMSolidSolver::createFromCube(double timeStep, double framePeriod)
{
	tetgenio cube;
	tetrahedralizeCube(cube);

	FEMSolidSolver* solver = new FEMSolidSolver(timeStep, framePeriod);

	for (int i = 0; i < cube.numberoftetrahedra; i++)
	{
		solver->tetraIndices.push_back(Eigen::Vector4i(cube.tetrahedronlist[4 * i + 0] - 1,
													   cube.tetrahedronlist[4 * i + 1] - 1,
													   cube.tetrahedronlist[4 * i + 2] - 1,
													   cube.tetrahedronlist[4 * i + 3] - 1));
	}

	for (int i = 0; i < cube.numberofpoints; i++)
	{
		Eigen::Vector3d zeroVec = Eigen::Vector3d(0.0, 0.0, 0.0);
		solver->positions.push_back(Eigen::Vector3d(cube.pointlist[3 * i + 0],
													cube.pointlist[3 * i + 1],
													cube.pointlist[3 * i + 2]));
		solver->velocities.push_back(zeroVec);

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
	Dm = std::vector<Eigen::Matrix3d>(tetraIndices.size());
	Bm = std::vector<Eigen::Matrix3d>(tetraIndices.size());
	We = std::vector<double>(tetraIndices.size());
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < tetraIndices.size(); i++)
	{
		Eigen::Matrix3d dmt;
		Eigen::Vector3d& lastCorner = positions[tetraIndices[i][3]];

		for (int cornerInd = 0; cornerInd < 3; ++cornerInd)
		{
			Eigen::Vector3d& corner = positions[tetraIndices[i][cornerInd]];
			dmt.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}

		Dm.at(i) = dmt;
		Bm.at(i) = dmt.inverse();
		We.at(i) = 1.0 / 6.0 * dmt.determinant();
	}
}

Eigen::Matrix3d FEMSolidSolver::computeP(Eigen::Matrix3d dst)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(dst, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV();
	Eigen::Vector3d Sigma = svd.singularValues();

	if (U.determinant() < 0)
	{
		U.col(2) *= -1.0;
	}
	if (V.determinant() < 0)
	{
		V.col(2) *= -1.0;
	}

	const double E = 0.1;
	const double NU = 0.4999;

	double mu = E / (2 * (1 + NU));
	double lambda = (E * NU) / ((1 + NU) * (1 - 2 * NU));

	Eigen::Matrix3d R = U*(V.transpose());
	Eigen::Matrix3d firstComponent = 2 * mu * (dst - R);

	Eigen::Matrix3d I = Eigen::MatrixXd::Identity(3, 3);
	Eigen::Matrix3d strainMeasurement = R * dst - I;
	Eigen::Matrix3d secondComponent = (lambda * (strainMeasurement.trace())) * R;

	Eigen::Matrix3d P = firstComponent + secondComponent;

	return P;

	/*const double E = 0.1;
	const double NU = 0.4999;

	double mu = E / (2 * (1 + NU));
	Eigen::Matrix3d P = 100.0 * (dst - Eigen::Matrix3d::Identity(3, 3));
	return P;*/
}

void FEMSolidSolver::computeBodyForce()
{
	for (size_t pointInd = 0; pointInd != positions.size(); ++pointInd)
	{
		bodyForces[pointInd] += Eigen::Vector3d(0.0, -masses.at(pointInd) * GravityAcc, 0.0);
	}
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

# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int pointInd = 0; pointInd < positions.size(); ++pointInd)
	{
		Eigen::Vector3d totalForce = bodyForces[pointInd] + elasticForces[pointInd];
		velocities[pointInd] += this->timeStep * (totalForce * (1.0 / masses[pointInd]));
		positions[pointInd] += this->timeStep * velocities[pointInd];

		if (positions[pointInd][1] < -1.0)
		{
			positions[pointInd][1] = -1.0;
			velocities[pointInd][1] = 0.0;
		}
	}

	if (this->steps % this->stepsPerFrame == 0)
	{
		std::string path = "./resultcache/polyCube";
		path += std::to_string(steps / stepsPerFrame);
		path += ".poly";
		this->save2File(path);
	}
	this->steps++;
}

void FEMSolidSolver::save2File(std::string path)
{
	std::fstream file;
	file.open(path, std::ios::out);
	objwriter::ObjWriter w(file);
	w.point();
	for (int i = 0; i < positions.size(); i++)
	{
		w.vertex(static_cast<float>(positions[i][0]),
				 static_cast<float>(positions[i][1]),
				 static_cast<float>(positions[i][2]), i);
	}
	w.polys();
	for (int i = 0; i < tetraIndices.size(); i++)
	{
		w.tet(tetraIndices[i][0] + 1, tetraIndices[i][1] + 1, tetraIndices[i][2] + 1, tetraIndices[i][3] + 1, i);
	}
	w.end();
	file.close();
}

void FEMSolidSolver::computeElasticForce()
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int tetInd = 0; tetInd < tetraIndices.size(); tetInd++)
	{
		Eigen::Matrix3d dst;
		Eigen::Vector3d& lastCorner = positions[tetraIndices[tetInd][3]];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Eigen::Vector3d& corner = positions[tetraIndices[tetInd][cornerInd]];
			dst.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}

		dst *= Bm[tetInd];
		
		Eigen::Matrix3d P = this->computeP(dst);

		// todo: what if rotation matrices are reflections
		if (dst.determinant() < 0)
		{

		}
		
		// populate H
		Eigen::Matrix3d H = -We[tetInd] * P * (Bm[tetInd].transpose());

		/// 注意！原版本存在初始化问题！
		Eigen::Vector3d f3 = Eigen::Vector3d(0.0, 0.0, 0.0);
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			elasticForces[tetraIndices[tetInd][cornerInd]] += H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		/// 等等……啥？应该是+=
		elasticForces[tetraIndices[tetInd][3]] += f3;
	}
}