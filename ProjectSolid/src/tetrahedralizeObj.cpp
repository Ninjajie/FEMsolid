//这个坑以后再填叭

# include "FEMSolid.h"
# include "tiny_obj_loader.h"

void tetrahedralizeObj(std::string objPath, tetgenio& out)
{
	tetgenio in;
	tetgenio::facet *f;
	tetgenio::polygon *p;

	in.firstnumber = 1;

	tinyobj::attrib_t attribs;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;

	bool loadStatus = tinyobj::LoadObj(&attribs, &shapes, &materials, &err, objPath.c_str());
	if (!err.empty())
	{
		std::cerr << err << std::endl;
	}

	tinyobj::shape_t& objShape = shapes[0];

	in.numberofpoints = attribs.vertices.size() / 3;
	in.pointlist = new REAL[in.numberofpoints * 3];
# ifdef OMParallelize
# pragma omp parallel for
# endif
	for (int pointInd = 0; pointInd < in.numberofpoints; ++pointInd)
	{
		in.pointlist[3 * pointInd + 0] = attribs.vertices.at(3 * pointInd + 0);
		in.pointlist[3 * pointInd + 1] = attribs.vertices.at(3 * pointInd + 1);
		in.pointlist[3 * pointInd + 2] = attribs.vertices.at(3 * pointInd + 2);
	}

	in.numberoffacets = objShape.mesh.num_face_vertices.size();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	size_t index_offset = 0;
	for (int faceInd = 0; faceInd < in.numberoffacets; ++faceInd)
	{
		f = &in.facetlist[faceInd];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = objShape.mesh.num_face_vertices[faceInd];
		p->vertexlist = new int[p->numberofvertices];
		for (int v = 0; v < p->numberofvertices; ++v)
		{
			tinyobj::index_t idx = objShape.mesh.indices[v + index_offset];
			p->vertexlist[v] = idx.vertex_index + 1;
		}
		index_offset += p->numberofvertices;
	}

	for (int facetInd = 0; facetInd < in.numberoffacets; ++facetInd)
	{
		in.facetmarkerlist[facetInd] = 0;
	}

	tetrahedralize("pq1.414a0.1", &in, &out);
}