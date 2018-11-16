#include "stdafx.h"
#include "Data_Structures.h"
#include "Data_Structures_Test.h"

#include <iostream>
//

bool Embedding_Dimension_Test_1()
{
	int num_points = 4;
	int dim = 2;
	double** p = new double*[num_points];
	for (int i = 0; i < num_points; i++)
		p[i] = new double[dim];

	p[0][0] = 0.0;
	p[0][1] = 0.0;
	p[1][0] = 1.0;
	p[1][1] = 0.0;
	p[2][0] = 0.0;
	p[2][1] = 1.0;
	p[3][0] = 1.0;
	p[3][1] = 1.0;

	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, p, num_points, dim);

	double** x = new double*[num_points];
	int Calculated_Dimension = Embedding_Dimension(Distances.Operator, num_points, x);
	//
	std::cout << "X : " << std::endl;
	for (int t1 = 0; t1 < num_points; t1++)
	{
		for (int t2 = 0; t2 < (num_points - 1); t2++)
		{
			std::cout << x[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//
	Distance_Struct Distances_New = Distance_Matrix(Euclidean_Distance, x, num_points, num_points - 1);

	std::cout << "Distance Matrix 0 : " << std::endl;
	std::cout << Distances.Operator << std::endl;
	std::cout << std::endl;
	std::cout << "Distance Matrix 1 : " << std::endl;
	std::cout << Distances_New.Operator << std::endl;
	std::cout << std::endl;

	return (Calculated_Dimension == dim);
}

bool Embedding_Dimension_Test_2()
{
	int num_points = 4;
	int dim = 2;
	double** p = new double*[num_points];
	for (int i = 0; i < num_points; i++)
		p[i] = new double[dim];

	p[0][0] = 1.35;
	p[0][1] = 2.65;
	p[1][0] = 1.44;
	p[1][1] = -2.97;
	p[2][0] = -0.53;
	p[2][1] = 4.28;
	p[3][0] = -0.87;
	p[3][1] = 4.76;

	//Distance Matrix
	//0        5.62072  2.48823  3.06276
	//5.62072        0  7.51288  8.06778
	//2.48823  7.51288        0  0.588218
	//3.06276  8.06778  0.588218        0
	//
	std::cout << "X0 : " << std::endl;
	for (int t1 = 0; t1 < num_points; t1++)
	{
		for (int t2 = 0; t2 < dim; t2++)
		{
			std::cout << p[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//

	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, p, num_points, dim);

	double** x = new double*[num_points];
	int Calculated_Dimension = Embedding_Dimension(Distances.Operator, num_points, x);
	//
	std::cout << "X : " << std::endl;
	for (int t1 = 0; t1 < num_points; t1++)
	{
		for (int t2 = 0; t2 < (num_points - 1); t2++)
		{
			std::cout << x[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//
	Distance_Struct Distances_New = Distance_Matrix(Euclidean_Distance, x, num_points, num_points - 1);

	std::cout << "Distance Matrix 0 : " << std::endl;
	std::cout << Distances.Operator << std::endl;
	std::cout << std::endl;
	std::cout << "Distance Matrix 1 : " << std::endl;
	std::cout << Distances_New.Operator << std::endl;
	std::cout << std::endl;

	return (Calculated_Dimension == dim);
}

bool Embedding_Dimension_Test_3()
{
	return Embedding_Dimension_Test_3(12);
}
bool Embedding_Dimension_Test_3(int num_points)
{
	int dim = 3;
	double** p = new double*[num_points];
	for (int i = 0; i < num_points; i++)
	{
		p[i] = new double[dim];
		for (int j = 0; j < dim; j++)
		{
			p[i][j] = 10.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 5.0;
		}
	}

	std::cout << "X0 : " << std::endl;
	for (int t1 = 0; t1 < num_points; t1++)
	{
		for (int t2 = 0; t2 < dim; t2++)
		{
			std::cout << p[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//

	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, p, num_points, dim);

	double** x = new double*[num_points];
	int Calculated_Dimension = Embedding_Dimension(Distances.Operator, num_points, x);
	//
	std::cout << "X : " << std::endl;
	for (int t1 = 0; t1 < num_points; t1++)
	{
		for (int t2 = 0; t2 < (num_points - 1); t2++)
		{
			std::cout << x[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//
	Distance_Struct Distances_New = Distance_Matrix(Euclidean_Distance, x, num_points, num_points - 1);

	std::cout << "Distance Matrix 0 : " << std::endl;
	std::cout << Distances.Operator << std::endl;
	std::cout << std::endl;
	std::cout << "Distance Matrix 1 : " << std::endl;
	std::cout << Distances_New.Operator << std::endl;
	std::cout << std::endl;

	return (Calculated_Dimension == dim);
}

bool Embedding_Dimension_Test_All()
{
	bool result_1 = Embedding_Dimension_Test_1();
	bool result_2 = Embedding_Dimension_Test_2();
	bool result_3 = Embedding_Dimension_Test_3();

	return (result_1 & result_2 & result_3);
}