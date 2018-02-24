# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

using namespace Eigen;
/*struct tetCorners
{
	std::array<Vector3d, 4> positions;
};*/

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
	int t;

	tetTetrahedras tetlist;
	tetCorners cornerlist;
	tetlist.tetNumber = t1.out.numberoftetrahedra;
	cornerlist.cornerNumber = t1.out.numberofpoints;

	for (int i = 0; i < tetlist.tetNumber; i++)
	{
		tetlist.indexes.push_back(Vector4i(t1.out.tetrahedronlist[4 * i],
										   t1.out.tetrahedronlist[4 * i + 1],
										   t1.out.tetrahedronlist[4 * i + 2],
										   t1.out.tetrahedronlist[4 * i + 3]));
	}

	for (int i = 0; i < cornerlist.cornerNumber; i++)
	{
		cornerlist.positions.push_back(Vector3d(t1.out.pointlist[3 * i],
												t1.out.pointlist[3 * i + 1],
												t1.out.pointlist[3 * i + 2]));
		cornerlist.velocities.push_back(Vector3d());
		cornerlist.elasticForces.push_back(Vector3d());
	}

	//=================================================================
	//Main Algorithm
	std::vector<Matrix3d> Dm;
	std::vector<Matrix3d> Bm;
	std::vector<double> We;

	//For all tet, precompute the volume
	for (int i = 0; i < tetlist.indexes.size(); i++)
	{
		Matrix3d dmt;
		Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[i][3]];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Vector3d& corner = cornerlist.positions[tetlist.indexes[i][cornerInd]];

			dmt.row(cornerInd) << corner - lastCorner;
		}
		Dm.push_back(dmt);
		Bm.push_back(dmt.inverse());
		We.push_back(dmt.determinant());
	}
	
	// For each frame, calculate new position
	int frameNumber = 1000;
	double deltaTime = 0.016;
	for (int i = 0; i < frameNumber; ++i)
	{
		// Renew velocity and position
		for (int j = 0; j < cornerlist.cornerNumber; i++)
		{
			cornerlist.velocities[j][1] -= 9.8 * deltaTime;
			cornerlist.positions[j] += cornerlist.velocities[i] * deltaTime;
			cornerlist.elasticForces[j] = Vector3d();
		}

		// Calculate elastic force
		for (int tetInd = 0; tetInd < t1.out.numberoftetrahedra; tetInd++)
		{
			Matrix3d dst;
			Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[tetInd][3]];
			for (int cornerInd = 0; cornerInd < 3; cornerInd++)
			{
				Vector3d& corner = cornerlist.positions[tetlist.indexes[tetInd][cornerInd]];
				dst.row(cornerInd) << corner - lastCorner;
			}
			dst *= Bm[tetInd];
			Matrix3d F;

			// todo: compute P
			Matrix3d p;

			// populate H
			Matrix3d H = -We[tetInd] * p * Bm[tetInd].transpose();

			Vector3d f3 = Vector3d();
			for (int cornerInd = 0; cornerInd < 3; cornerInd++)
			{
				cornerlist.elasticForces[tetlist.indexes[tetInd][cornerInd]] += H.col(cornerInd);
				f3 -= H.col(cornerInd);
			}
			cornerlist.elasticForces[tetlist.indexes[tetInd][3]] += f3;
		}

		// Boundary situation
		for (int j = 0; j < cornerlist.cornerNumber; i++)
		{
			cornerlist.velocities[j] += cornerlist.elasticForces[i];
			if (cornerlist.positions[j][1] < -1.)
			{
				cornerlist.positions[j][1] = -1.;
				cornerlist.velocities[j][1] = 0.8 * fabs(cornerlist.velocities[i][1]);
			}
		}
	}
}