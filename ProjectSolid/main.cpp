# include <iostream>
# include "FEMSolid.h"
# include "tet.h"

using namespace Eigen;
struct tetH
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
	std::vector<tetH> tetlst;
	int  t;
	std::cout<<t1.out.numberoftetrahedra<<std::endl;
	std::cout << t1.out.numberofcorners<<std::endl;
	for (int i = 0; i < t1.out.numberoftetrahedra; i++)
	{
		tetH temp;
		
		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 3; k++)
			{
				temp.positions[j][k]=t1.out.pointlist[3 * (t1.out.tetrahedronlist[4 * i + j]) + k];
			}

		}
		tetlst.push_back(temp);
	}

	for (int i = 0; i < tetlst.size(); i++)
	{
		std::cout << "the" << i << "th tetrahedral"<<std::endl;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				std::cout << tetlst[i].positions[j][k];
				std::cout << "||";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cin >> t;
	//=================================================================
}