# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

#define E 0.1
#define NU 0.4999
using namespace Eigen;

struct tetTetrahedras
{
	int tetNumber;
	std::vector<Vector4i> indexes;
};

struct tetCorners
{
	int cornerNumber;
	std::vector<Vector3d> positions;
	std::vector<Vector3d> velocities;
	std::vector<Vector3d> elasticForces;
};

int main(int argc, char** argv)
{
	std::cout << "Just seeing if this works" << std::endl;
	Eigen::Matrix4f mat = Eigen::Matrix4f();
	std::cout << mat << std::endl;
	tetgenio ioOp = tetgenio();
	ioOp.initialize();
	//=========================================== test tetGen
	
	tetGen t1;

	tetTetrahedras tetlist;
	tetCorners cornerlist;
	tetlist.tetNumber = t1.out.numberoftetrahedra;
	cornerlist.cornerNumber = t1.out.numberofpoints;

	for (int i = 0; i < tetlist.tetNumber; i++)
	{
		tetlist.indexes.push_back(Vector4i(t1.out.tetrahedronlist[4 * i] - 1,
										   t1.out.tetrahedronlist[4 * i + 1] - 1,
										   t1.out.tetrahedronlist[4 * i + 2] - 1,
										   t1.out.tetrahedronlist[4 * i + 3] - 1));
	}

	for (int i = 0; i < cornerlist.cornerNumber; i++)
	{
		cornerlist.positions.push_back(Vector3d(t1.out.pointlist[3 * i],
												t1.out.pointlist[3 * i + 1],
												t1.out.pointlist[3 * i + 2]));
		cornerlist.velocities.push_back(Vector3d(0, 0, 0));
		cornerlist.elasticForces.push_back(Vector3d(0, 0, 0));
	}

	//=================================================================
	//Main Algorithm
	std::vector<Matrix3d> Dm;
	std::vector<Matrix3d> Bm;
	std::vector<double> We;

	//For all tet, precompute the volume
	for (int i = 0; i < tetlist.tetNumber; i++)
	{
		Matrix3d dmt;
		Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[i][3]];
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Vector3d& corner = cornerlist.positions[tetlist.indexes[i][cornerInd]];

			dmt.col(cornerInd) << corner[0] - lastCorner[0] , corner[1] - lastCorner[1], corner[2] - lastCorner[2];
		}
		Dm.push_back(dmt);
		Bm.push_back(dmt.inverse());
		We.push_back(dmt.determinant());
	}
	
	// For each frame, calculate new position
	int frameNumber = 30;
	double deltaTime = 0.016;
	for (int i = 0; i < frameNumber; ++i)
	{
		// Renew velocity and position
		for (int j = 0; j < cornerlist.cornerNumber; ++j)
		{
			cornerlist.velocities[j][1] -= 9.8 * deltaTime;

			//std::cout << cornerlist.velocities[j] << std::endl;
			cornerlist.positions[j] += cornerlist.velocities[j] * deltaTime;
			std::cout << cornerlist.positions[j] << std::endl;
			cornerlist.elasticForces[j] = Vector3d(0, 0, 0);
		}

		// Calculate elastic force
		for (int tetInd = 0; tetInd < tetlist.tetNumber; tetInd++)
		{
			Matrix3d dst;
			Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[tetInd][3]];
			for (int cornerInd = 0; cornerInd < 3; cornerInd++)
			{
				Vector3d& corner = cornerlist.positions[tetlist.indexes[tetInd][cornerInd]];
				dst.col(cornerInd) << corner[0] - lastCorner[0], corner[1] - lastCorner[1], corner[2] - lastCorner[2];
			}
			dst *= Bm[tetInd];

			Matrix3d p;
			JacobiSVD<Matrix3d> svd(dst, ComputeFullU | ComputeFullV);
			Matrix3d U = svd.matrixU();
			Matrix3d V = svd.matrixV();
			Vector3d Sigma = svd.singularValues();
			double mu = E / (2 * (1 + NU));
			double lambda = (E * NU) / ((1 + NU) * (1 - 2 * NU));
			Matrix3d R = U*V;
			p = 2 * mu * (dst - R) + lambda * (R.transpose() * dst - MatrixXd::Identity(3, 3)).trace() * R;
			// todo: what if rotation matrices are reflections
			if (dst.determinant() < 0)
			{

			}

			// populate H
			Matrix3d H = -We[tetInd] * p * Bm[tetInd].transpose();

			Vector3d f3 = Vector3d(0, 0, 0);
			for (int cornerInd = 0; cornerInd < 3; cornerInd++)
			{
				cornerlist.elasticForces[tetlist.indexes[tetInd][cornerInd]] += H.col(cornerInd);
				f3 -= H.col(cornerInd);
			}
			cornerlist.elasticForces[tetlist.indexes[tetInd][3]] += f3;
		}

		// Boundary situation
		for (int j = 0; j < cornerlist.cornerNumber; ++j)
		{
			cornerlist.velocities[j] += cornerlist.elasticForces[j] * deltaTime * 0.001;
			if (cornerlist.positions[j][1] < -1.)
			{
				cornerlist.positions[j][1] = -1.;
				cornerlist.velocities[j][1] = 0.8 * fabs(cornerlist.velocities[j][1]);
			}
		}


		//=======================================test output===========

		std::fstream file;
		file.open("test" + std::to_string(i) + ".poly", std::ios::out);
		objwriter::ObjWriter w(file);
		w.point();
		for (int i = 0; i < cornerlist.positions.size(); i++)
		{
			w.vertex(float(cornerlist.positions[i][0]), float(cornerlist.positions[i][1]), float(cornerlist.positions[i][2]), i);
		}
		w.polys();
		for (int i = 0; i < tetlist.indexes.size(); i++)
		{
			w.tet(tetlist.indexes[i][0], tetlist.indexes[i][1], tetlist.indexes[i][2], tetlist.indexes[i][3], i);
		}
		w.end();
		//=============================================================
	}
}