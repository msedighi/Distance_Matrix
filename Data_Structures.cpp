#include "stdafx.h"
#include "Data_Structures.h"

//
double Euclidean_Distance(double* v1, double* v2, int dim)
{
	double norm = 0;
	for (int i = 0; i < dim; i++)
		norm += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	return sqrt(norm);
}

double Euclidean_Distance_Vector(Eigen::VectorXd v1, Eigen::VectorXd v2)
{
	Eigen::VectorXd dv = v1 - v2;
	return dv.norm();
}

template<int p>
double lpDistance(Eigen::VectorXd v1, Eigen::VectorXd v2)
{
	return (v1 - v2).lpNorm<p>();
}

int compareDistance(const void * a, const void * b)
{
	if ((*(Distance_Index*)a).Distance < (*(Distance_Index*)b).Distance) return -1;
	if ((*(Distance_Index*)a).Distance == (*(Distance_Index*)b).Distance) return 0;
	if ((*(Distance_Index*)a).Distance > (*(Distance_Index*)b).Distance) return 1;
}

template<class MyType>
int compareMyType(const void * a, const void * b)
{
	if (*(MyType*)a < *(MyType*)b) return -1;
	if (*(MyType*)a == *(MyType*)b) return 0;
	if (*(MyType*)a > *(MyType*)b) return 1;
}
//

double*** Difference_Matrix(double** vectors, int num_points, int dim)
{
	double*** difference_array = new double**[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		difference_array[index_r] = new double*[num_points];
		for (int index_c = 0; index_c < num_points; index_c++)
		{
			for (int i = 0; i < dim; i++) difference_array[index_r][index_c][i] = vectors[index_c][i] - vectors[index_r][i];
		}
	}
	return difference_array;
}

double* Vector_Product(double** x, double** y, int num_points, int dim)
{
	double* product = new double[num_points];

	for (int index_r = 0; index_r < num_points; index_r++)
	{
		product[index_r] = 0;
		for (int index_c = 0; index_c < dim; index_c++)
		{
			product[index_r] += x[index_r][index_c] * y[index_r][index_c];
		}
	}

	return product;
}

Distance_Struct Distance_Matrix(double(*Distance_func)(double*, double*, int), double** Positions, long Number_Points, int dim)
{
	long Number_Pairs = Number_Points * (Number_Points - 1) / 2;
	Distance_Struct distance = Distance_Struct(Number_Points);

	int counter = 0;
	for (long index_r = 0; index_r < Number_Points; index_r++)
	{
		for (long index_c = index_r + 1; index_c < Number_Points; index_c++)
		{
			distance.Operator(index_r, index_c) = Distance_func(Positions[index_r], Positions[index_c], dim);
			distance.Operator(index_c, index_r) = distance.Operator(index_r, index_c);
			distance.Vector[counter].Index1 = index_r;
			distance.Vector[counter].Index2 = index_c;
			distance.Vector[counter].Distance = distance.Operator(index_r, index_c);
			counter++;
		}
	}

	qsort(distance.Vector, Number_Pairs, sizeof(Distance_Index), compareDistance);

	return distance;
}

int Embedding_Dimension(Eigen::MatrixXd distance_matrix, int num_points, double** x)
{
	x = new double*[num_points];

	int dim = 1;
	x[0] = new double[1];
	x[1] = new double[1];
	x[0][0] = 0.0;
	x[1][0] = distance_matrix(0, 1);
	for (int j_p = 2; j_p < num_points; j_p++)
	{
		for (int i_p = 0; i_p < j_p; i_p++)
		{
			double* a = Solver(x, dim + 1, distance_matrix.row(dim + 1));
			if (isnan(a[dim]))
				return INFINITY;
			else
			{
				for (int i = 0; i < dim; i++)
				{
					x[i][dim] = 0.0;
					x[dim + 1][i] = a[i];
				}
				dim++;
			}
		}
	}

	return dim;
}

double* Solver(double** x, int dim, Eigen::VectorXd distance_vector)
{

	// "distance_vector" has (dim) elements containing the distance between the (dim + 1) point & the other (dim) points 
	// "x" is (dim, dim) giving the coordinates of the (dim) points in the current dimension, (dim)
	// So we have (dim) equations & (dim) unknown coordinates
	double norm2_x_new = distance_vector[0] * distance_vector[0];
	double* a = new double[dim];
	for (int i = 0; i < (dim - 1); i++)
	{
		double norm2_x = 0;
		for (int i_d = 0; i_d < (dim - 1); i_d++)
		{
			norm2_x += x[i + 1][i_d] * x[i + 1][i_d];
		}

		a[i] = (norm2_x + norm2_x_new - distance_vector[i + 1] * distance_vector[i + 1]) / 2.0;
	}

	for (int i = 1; i < dim; i++)
	{
		if (abs(x[i][i - 1]) > 0)
			a[i-1] = a[i-1] / x[i][i - 1];
		else
			break;

		for (int j = (i + 1); j < dim; j++)
		{
			a[j-1] = (a[j-1] - x[j][i - 1] * a[i-1]) / x[i][i - 1];
		}
	}

	double norm2_a = 0;
	for (int i = 0; i < (dim - 1); i++)
	{
		norm2_a += a[i] * a[i];
	}
	if (norm2_x_new > norm2_a)
	{
		a[dim - 1] = sqrt(norm2_x_new - norm2_a);

	}
	else
		a[dim - 1] = NAN;

	return a;
}