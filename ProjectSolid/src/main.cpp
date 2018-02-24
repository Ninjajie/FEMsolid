# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

#define E 0.1
#define NU 0.4999
using namespace Eigen;
/*struct tetCorners
{
	std::array<Vector3d, 4> positions;
};*/

struct tetTetrahedras
{
	std::vector<Vector4i> indexes;
};

struct tetCorners
{
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
	int t;
	std::cout << t1.out.numberoftetrahedra << std::endl;
	std::cout << t1.out.numberofcorners << std::endl;

	for (int i = 0; i < t1.out.numberoftetrahedra; i++)
	{
		tetlist.indexes.push_back(Vector4i(t1.out.tetrahedronlist[4 * i],
										   t1.out.tetrahedronlist[4 * i + 1],
										   t1.out.tetrahedronlist[4 * i + 2],
										   t1.out.tetrahedronlist[4 * i + 3]));
	}

	for (int i = 0; i < t1.out.numberofcorners; i++)
	{
		Vector3d zero = Vector3d();
		cornerlist.positions.push_back(Vector3d(t1.out.pointlist[3 * i],
												t1.out.pointlist[3 * i + 1],
												t1.out.pointlist[3 * i + 2]));
		cornerlist.velocities.push_back(zero);
		cornerlist.elasticForces.push_back(zero);
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
	// Compute the elastic forces
	// current positions
	// todo: need to populate it
	std::vector<tetCorners> currentTetList;
	
	// force on each corner
	std::vector<tetCorners> tetForces;
	//stress tensor
	//std::vector<Matrix3d> F;

	for (int tetInd = 0; tetInd < tetlist.indexes.size(); tetInd++)
	{
		Matrix3d dst;
		Vector3d& lastCorner = cornerlist.positions[tetlist.indexes[tetInd][3]];
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Vector3d& corner = cornerlist.positions[tetlist.indexes[tetInd][cornerInd]];

			dst.row(cornerInd) << corner - lastCorner;
		}
		dst *= Bm[tetInd];
		// todo: compute P
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
		Matrix3d H = - We[tetInd] * p * Bm[tetInd].transpose();

		Vector3d f3 = Vector3d();
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			cornerlist.elasticForces[tetlist.indexes[tetInd][cornerInd]] = H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		
		cornerlist.elasticForces[tetlist.indexes[tetInd][3]] = f3;
	}
}