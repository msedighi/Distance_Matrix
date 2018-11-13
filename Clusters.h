#pragma once
#ifndef CLUSTERS_H
#define CLUSTERS_H
#endif // !1

#include "Data_Structures.h"

class Clusters
{
private:
	void Hierarchical(Distance_Struct Distances);
public:
	Eigen::MatrixXd Orthogonal_Transformation;
	Eigen::MatrixXd Orthogonal_Transformation_Dual;
	Dendo_Slice * Dendogram;
	Dendo_Slice * Dendogram_Dual;
	double Max_Vacuum_Scales;
	double Min_Saturation_Scales;
		
	Clusters(Distance_Struct Distances);
	~Clusters();
};

