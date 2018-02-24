#define TETLIBRARY
#ifndef TET_H
#define TET_H
#include"FEMSolid.h"

class tetGen
{
public:
	tetGen();
	~tetGen() {}
	tetgenio in, out;
	tetgenio::facet *f;
	tetgenio::polygon *p;

};

#endif // !TET_H
