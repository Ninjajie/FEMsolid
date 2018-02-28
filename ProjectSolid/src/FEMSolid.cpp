# include "FEMSolid.h"

void tetrahedralizeCube(tetgenio& out);
void tetrahedralizeObj(std::string objPath, tetgenio& out);

FEMSolidSolver::FEMSolidSolver(tetgenio& mesh, fReal timeStep, fReal framePeriod)
	: timeStep(timeStep), framePeriod(framePeriod), steps(0), stepsPerFrame(static_cast<long long>(framePeriod / timeStep)), frames(0)
{
# ifdef OMParallelize
	omp_set_num_threads(TOTALThreads);
# endif
	this->preAllocate(mesh.numberoftetrahedra, mesh.numberofpoints);
	this->preFill(mesh);
	this->setInitial();
	this->preCompute();
}

long long FEMSolidSolver::getCurrentIterations()
{
	return this->steps;
}

long long FEMSolidSolver::getCurrentFrames()
{
	return this->frames;
}

void FEMSolidSolver::preAllocate(int numOfTets, int numOfVerts)
{
	this->numOfTets = numOfTets;
	this->numOfVerts = numOfVerts;

	this->tetraIndices = new Eigen::Vector4i[numOfTets];
	this->Dm = new mat3[numOfTets];
	this->Bm = new mat3[numOfTets];
	this->We = new fReal[numOfTets];
	
	this->positions = new vec3[numOfVerts];
	this->velocities = new vec3[numOfVerts];
	this->elasticForces = new vec3[numOfVerts];
	
	this->bodyForces = new vec3[numOfVerts];
	this->masses = new fReal[numOfVerts];
}

FEMSolidSolver::~FEMSolidSolver()
{
	delete[] this->tetraIndices;
	delete[] this->Dm;
	delete[] this->Bm;
	delete[] this->We;

	delete[] this->positions;
	delete[] this->velocities;
	delete[] this->elasticForces;

	delete[] this->bodyForces;
	delete[] this->masses;
}

void FEMSolidSolver::preFill(tetgenio& mesh)
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < numOfTets; i++)
	{
		tetraIndices[i] = Eigen::Vector4i(mesh.tetrahedronlist[4 * i + 0] - 1,
										  mesh.tetrahedronlist[4 * i + 1] - 1,
										  mesh.tetrahedronlist[4 * i + 2] - 1,
										  mesh.tetrahedronlist[4 * i + 3] - 1);
	}
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < numOfVerts; i++)
	{
		positions[i] = vec3(mesh.pointlist[3 * i + 0],
							mesh.pointlist[3 * i + 1],
							mesh.pointlist[3 * i + 2]);
		
	}
}

void FEMSolidSolver::setInitial()
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < numOfVerts; i++)
	{
		velocities[i] = zeroVec;
	}
}

FEMSolidSolver* FEMSolidSolver::createFromCube(fReal timeStep, fReal framePeriod)
{
	tetgenio cube;
	tetrahedralizeCube(cube);
	FEMSolidSolver* solver = new FEMSolidSolver(cube, timeStep, framePeriod);
	
	return solver;
}

FEMSolidSolver* FEMSolidSolver::createFromObj(const std::string& filename, fReal timeStep, fReal framePeriod)
{
	tetgenio obj;
	tetrahedralizeObj(filename, obj);
	FEMSolidSolver* solver = new FEMSolidSolver(obj, timeStep, framePeriod);

	return solver;
}

void FEMSolidSolver::preCompute()
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < this->numOfVerts; ++i)
	{
		masses[i] = 0.0;
	}
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int i = 0; i < this->numOfTets; i++)
	{
		mat3 dmt;
		vec3& lastCorner = positions[tetraIndices[i][3]];

		for (int cornerInd = 0; cornerInd < 3; ++cornerInd)
		{
			vec3 corner = positions[tetraIndices[i][cornerInd]];
			dmt.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}

		Dm[i] = dmt;
		Bm[i] = dmt.inverse();
		fReal volume = 1.0 / 6.0 * std::abs(dmt.determinant());
		We[i] = volume;
		fReal tetMass = volume * Density;
		for (int vert = 0; vert != 4; ++vert)
		{
			masses[tetraIndices[i][vert]] += 0.25 * tetMass;
		}
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

	fReal mu = E / (2 * (1 + NU));
	fReal lambda = (E * NU) / ((1 + NU) * (1 - 2 * NU));

	mat3 R = U*(V.transpose());
	mat3 firstComponent = 2 * mu * (dst - R);

	mat3 I = mat3::Identity(3, 3);
	mat3 strainMeasurement = R.transpose() * dst - I;
	mat3 secondComponent = (lambda * (strainMeasurement.trace())) * R;

	mat3 P = firstComponent + secondComponent;

	return P;
}

void FEMSolidSolver::computeBodyForce()
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int pointInd = 0; pointInd < this->numOfVerts; ++pointInd)
	{
		bodyForces[pointInd] += vec3(0.0, -masses[pointInd] * GravityAcc, 0.0);
	}
}

void FEMSolidSolver::stepForward()
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int pointInd = 0; pointInd < this->numOfVerts; ++pointInd)
	{
		bodyForces[pointInd] = (zeroVec);
		elasticForces[pointInd] = (zeroVec);
	}

	this->computeBodyForce();
	this->computeElasticForce();

# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int pointInd = 0; pointInd < this->numOfVerts; ++pointInd)
	{
		vec3 totalForce = bodyForces[pointInd] + elasticForces[pointInd];
		velocities[pointInd] += this->timeStep * (totalForce * (1.0 / masses[pointInd]));
		solveForBoundary(pointInd);
		positions[pointInd] += this->timeStep * velocities[pointInd];
	}

	if (this->steps % this->stepsPerFrame == 0)
	{
		std::string path = "./resultcache/polyCube";
		path += std::to_string(steps / stepsPerFrame);
		path += ".poly";
		this->save2File(path);
		this->frames++;
	}
	this->steps++;
}

void FEMSolidSolver::solveForBoundary(int pointInd)
{
	if (std::abs(positions[pointInd][2]) < 1e-3)
	{
		positions[pointInd][2] = 0.0;
		velocities[pointInd] = zeroVec;
	}
	/*if (positions[pointInd][1] < -1.0)
	{
		positions[pointInd][1] = -1.0;
		velocities[pointInd][1] = 0.0;
	}*/
}

void FEMSolidSolver::save2File(std::string path)
{
	std::fstream file;
	file.open(path, std::ios::out);
	objwriter::ObjWriter w(file);
	w.point();
	for (int i = 0; i < this->numOfVerts; i++)
	{
		w.vertex(static_cast<fReal>(positions[i][0]),
				 static_cast<fReal>(positions[i][1]),
				 static_cast<fReal>(positions[i][2]), i);
	}
	w.polys();
	for (int i = 0; i < this->numOfTets; i++)
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
	for (int tetInd = 0; tetInd < this->numOfTets; tetInd++)
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
		vec3 f3 = zeroVec;
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
# ifdef OMParallelize
# pragma omp critical
# endif
			elasticForces[tetraIndices[tetInd][cornerInd]] = elasticForces[tetraIndices[tetInd][cornerInd]] + H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		/// 等等……啥？应该是+=
# ifdef OMParallelize
# pragma omp critical
# endif
		elasticForces[tetraIndices[tetInd][3]] = elasticForces[tetraIndices[tetInd][3]] + f3;
	}
}