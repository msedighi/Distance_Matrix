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
	Dendo * Dendogram;
	Dendo * Dendogram_Dual;
	double Max_Vacuum_Scale;
	double Min_Saturation_Scale;

	Clusters(Distance_Struct Distances);
	~Clusters();
};

