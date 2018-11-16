#include "stdafx.h"
#include "Eigen/Dense"
#include "Clusters.h"

#include <iostream>

using namespace Eigen;

Clusters::Clusters(Distance_Struct Distances)
{
	Orthogonal_Transformation = MatrixXd::Constant(Distances.Number_Points, Distances.Number_Points, 1 / sqrt(Distances.Number_Points)); // Output
	Orthogonal_Transformation_Dual = MatrixXd::Constant(Distances.Number_Points, Distances.Number_Points, 1 / sqrt(Distances.Number_Points)); // Output

	// Defining the Boolean Dendogram
	Dendogram = new Dendo(Distances.Number_Points);
	Dendogram_Dual = new Dendo(Distances.Number_Points);

	// Initializing the base of the Dendogram
	Dendogram->Scale[0] = 0;
	Dendogram_Dual->Scale[0] = Infinity;

	Hierarchical(Distances);

	Max_Vacuum_Scale = Dendogram->Scale[Distances.Number_Points - 1];
	Min_Saturation_Scale = Dendogram_Dual->Scale[Distances.Number_Points - 1];
}

void Clusters::Hierarchical(Distance_Struct Distances)
{
	// NOTE!!!! We need to also add the option for Hierarchical Clusterings from a local POV! *****

	// Defining the a "Cluster State" for each scale 
	ArrayXb * Clusters_State = new ArrayXb[Distances.Number_Points];
	int * Cluster_ID = new int[Distances.Number_Points];

	double Num1, Num2;

	// The Original Model (Small to Large Scale) 

	// Initiatation
	for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
	{
		Clusters_State[i_points] = ArrayXb::Constant(Distances.Number_Points, false);
		Clusters_State[i_points](i_points) = true;

		Cluster_ID[i_points] = -i_points;
	}
	unsigned cluster_counter = 0;
	//
	for (int i_pairs = 0; i_pairs < Distances.Number_Pairs; i_pairs++)
	{
		if (Cluster_ID[Distances.Vector[i_pairs].Index1] != Cluster_ID[Distances.Vector[i_pairs].Index2])
		{
			cluster_counter++;
			Num1 = Clusters_State[Distances.Vector[i_pairs].Index1].count();
			Num2 = Clusters_State[Distances.Vector[i_pairs].Index2].count();
			Orthogonal_Transformation.row(cluster_counter - 1) =
				(Clusters_State[Distances.Vector[i_pairs].Index1].cast<double>() * sqrt(Num2 / Num1)
					- Clusters_State[Distances.Vector[i_pairs].Index2].cast<double>() * sqrt(Num1 / Num2))
				/ sqrt(Num1 + Num2);

			// *** Dendogram Assignment
			ArrayXb this_Cluster_State = Clusters_State[Distances.Vector[i_pairs].Index1] || Clusters_State[Distances.Vector[i_pairs].Index2];
			
			for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
			{				 
				Dendogram->Structure[cluster_counter][i_points] = this_Cluster_State(i_points);
			}
			Dendogram->Scale[cluster_counter] = Distances.Vector[i_pairs].Distance;
			// ***
			for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
			{
				if (this_Cluster_State(i_points))
				{
					Cluster_ID[i_points] = cluster_counter;
					Clusters_State[i_points] = this_Cluster_State;
				}
			}

		}
	}
	// The Dual Model (Large to Small Scale) 

	// Initiatation
	for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
	{
		Clusters_State[i_points] = ArrayXb::Constant(Distances.Number_Points, false);
		Clusters_State[i_points](i_points) = true;

		Cluster_ID[i_points] = -i_points;
	}
	cluster_counter = 0;
	//
	for (int i_pairs = (Distances.Number_Pairs - 1); i_pairs >= 0; i_pairs--)
	{
		if (Cluster_ID[Distances.Vector[i_pairs].Index1] != Cluster_ID[Distances.Vector[i_pairs].Index2])
		{
			cluster_counter++;
			Num1 = Clusters_State[Distances.Vector[i_pairs].Index1].count();
			Num2 = Clusters_State[Distances.Vector[i_pairs].Index2].count();
			Orthogonal_Transformation_Dual.row(cluster_counter - 1) =
				(Clusters_State[Distances.Vector[i_pairs].Index1].cast<double>() * sqrt(Num2 / Num1)
					- Clusters_State[Distances.Vector[i_pairs].Index2].cast<double>() * sqrt(Num1 / Num2))
				/ sqrt(Num1 + Num2);

			// *** Dendogram Assignment
			ArrayXb this_Cluster_State = Clusters_State[Distances.Vector[i_pairs].Index1] || Clusters_State[Distances.Vector[i_pairs].Index2];
			for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
			{
				Dendogram_Dual->Structure[cluster_counter][i_points] = this_Cluster_State(i_points);
			}
			Dendogram_Dual->Scale[cluster_counter] = Distances.Vector[i_pairs].Distance;
			// ***
			for (int i_points = 0; i_points < Distances.Number_Points; i_points++)
			{
				if (this_Cluster_State(i_points))
				{
					Cluster_ID[i_points] = cluster_counter;
					Clusters_State[i_points] = this_Cluster_State;
				}
			}

		}
	}
}

Clusters::~Clusters() 
{
	delete Dendogram;
	delete Dendogram_Dual;
}
