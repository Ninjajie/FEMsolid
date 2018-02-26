# include "FEMSolid.h"

void tetrahedralizeCube(tetgenio& out);

FEMSolidSolver::FEMSolidSolver(fReal timeStep, fReal framePeriod)
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

FEMSolidSolver* FEMSolidSolver::createFromCube(fReal timeStep, fReal framePeriod)
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
		vec3 zeroVec = vec3(0.0, 0.0, 0.0);
		solver->positions.push_back(vec3(cube.pointlist[3 * i + 0],
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
	Dm = std::vector<mat3>(tetraIndices.size());
	Bm = std::vector<mat3>(tetraIndices.size());
	We = std::vector<fReal>(tetraIndices.size());
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < tetraIndices.size(); i++)
	{
		mat3 dmt;
		vec3& lastCorner = positions[tetraIndices[i][3]];

		for (int cornerInd = 0; cornerInd < 3; ++cornerInd)
		{
			vec3 corner = positions[tetraIndices[i][cornerInd]];
			dmt.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}

		Dm.at(i) = dmt;
		Bm.at(i) = dmt.inverse();
		We.at(i) = 1.0 / 6.0 * std::abs(dmt.determinant());
	}
}

mat3 FEMSolidSolver::computeP(mat3 dst)
{
	Eigen::JacobiSVD<mat3> svd(dst, Eigen::ComputeFullU | Eigen::ComputeFullV);
	mat3 U = svd.matrixU();
	mat3 V = svd.matrixV();
	//mat3 Sigma = svd.singularValues();

	if (U.determinant() < 0)
	{
		U.col(2) *= -1.0;
	}
	if (V.determinant() < 0)
	{
		V.col(2) *= -1.0;
	}

	const fReal E = 1000;
	const fReal NU = 0.3;

	fReal mu = E / (2 * (1 + NU));
	fReal lambda = (E * NU) / ((1 + NU) * (1 - 2 * NU));

	mat3 R = U*(V.transpose());
	mat3 firstComponent = 2 * mu * (dst - R);

	mat3 I = mat3::Identity(3, 3);
	mat3 strainMeasurement = R * dst - I;
	mat3 secondComponent = (lambda * (strainMeasurement.trace())) * R;

	mat3 P = firstComponent + secondComponent;

	return P;

	/*const fReal E = 0.1;
	const fReal NU = 0.4999;

	fReal mu = E / (2 * (1 + NU));
	mat3 P = 100.0 * (dst - mat3::Identity(3, 3));
	return P;*/
}

void FEMSolidSolver::computeBodyForce()
{
	for (size_t pointInd = 0; pointInd != positions.size(); ++pointInd)
	{
		bodyForces[pointInd] += vec3(0.0, -masses.at(pointInd) * GravityAcc, 0.0);
	}
}

void FEMSolidSolver::stepForward()
{
	vec3 zeroVec = vec3(0.0, 0.0, 0.0);
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
		vec3 totalForce = bodyForces[pointInd] + elasticForces[pointInd];
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
		w.vertex(static_cast<fReal>(positions[i][0]),
				 static_cast<fReal>(positions[i][1]),
				 static_cast<fReal>(positions[i][2]), i);
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
		mat3 dst;
		vec3& lastCorner = positions[tetraIndices[tetInd][3]];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			vec3& corner = positions[tetraIndices[tetInd][cornerInd]];
			dst.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}

		dst *= Bm[tetInd];
		
		mat3 P = this->computeP(dst);

		// todo: what if rotation matrices are reflections
		if (dst.determinant() < 0)
		{

		}
		
		// populate H
		mat3 H = -We[tetInd] * P * (Bm[tetInd].transpose());

		/// 注意！原版本存在初始化问题！
		vec3 f3 = vec3(0.0, 0.0, 0.0);
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			elasticForces[tetraIndices[tetInd][cornerInd]] += H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		/// 等等……啥？应该是+=
		elasticForces[tetraIndices[tetInd][3]] += f3;
	}
}