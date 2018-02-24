# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

using namespace Eigen;
struct tetCorners
{
	std::array<Vector3d, 4> positions;

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
	std::vector<tetCorners> tetlist;
	int t;
	std::cout << t1.out.numberoftetrahedra << std::endl;
	std::cout << t1.out.numberofcorners << std::endl;

	for (int i = 0; i < t1.out.numberoftetrahedra; i++)
	{
		tetCorners temp;
		
		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 3; k++)
			{
				temp.positions[j][k]=t1.out.pointlist[3 * (t1.out.tetrahedronlist[4 * i + j]) + k];
			}

		}
		tetlist.push_back(temp);
	}
	/*
	for (int i = 0; i < tetlist.size(); i++)
	{
		std::cout << " the " << i << "th tetrahedral"<<std::endl;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				std::cout << tetlist[i].positions[j][k];
				std::cout << " || ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cin >> t;
	*/
	//=================================================================
	//Main Algorithm
	std::vector<Matrix3d> Dm;
	std::vector<Matrix3d> Bm;
	std::vector<double> We;

	//For all tet, precompute the volume
	for (int i = 0; i < tetlist.size(); i++)
	{
		Matrix3d dmt;
		Vector3d& lastCorner = tetlist[i].positions[3];

		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Vector3d& corner = tetlist[i].positions[cornerInd];

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
	std::vector<Matrix3d> F;
	//std::vector<Vector3d> H;
	//std::vector<Matrix3d> P;
	for (int tetInd = 0; tetInd < tetlist.size(); tetInd++)
	{
		Matrix3d dst;
		Vector3d& lastCorner = currentTetList[tetInd].positions[3];
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			Vector3d& corner = currentTetList[tetInd].positions[cornerInd];

			dst.row(cornerInd) << corner - lastCorner;
		}
		dst *= Bm[tetInd];
		F.push_back(dst);

		// todo: compute P
		JacobiSVD<Matrix3d> svd(dst, ComputeFullU | ComputeFullV);
		Matrix3d p;

		// populate H
		Matrix3d H = - We[tetInd] * p * Bm[tetInd].transpose();
		tetCorners forcesFromCurrTet;

		Vector3d f3;
		f3 << 0, 0, 0;
		for (int cornerInd = 0; cornerInd < 3; cornerInd++)
		{
			forcesFromCurrTet.positions[cornerInd] = H.col(cornerInd);
			f3 -= H.col(cornerInd);
		}
		
		forcesFromCurrTet.positions[3] = f3;	
		tetForces.push_back(forcesFromCurrTet);
	}
}