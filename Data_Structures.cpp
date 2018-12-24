#include "stdafx.h"
#include "Data_Structures.h"
//
Dendo::Dendo(int num_points)
{
	Num_Points = num_points;

	Structure = new bool*[num_points];
	Scale = new double[num_points];
	for (int i = 0; i < num_points; i++)
	{
		Structure[i] = new bool[num_points];
		for (int j = 0; j < num_points; j++)
			Structure[i][j] = false;
	} 
}

Dendo::~Dendo()
{
	delete[] Scale;
	for (int i = 0; i < Num_Points; i++)
	{
		delete[] Structure[i];
	}
	delete [] Structure;
}

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
	for (int i = 0; i < num_points; i++)
	{
		x[i] = new double[num_points - 1];
		for (int j = 0; j < (num_points - 1); j++)
			x[i][j] = 0.0;
	}

	int dim = 1;
	x[1][0] = distance_matrix(0, 1);
	for (int i_p = 2; i_p < num_points; i_p++)
	{
		double* a = Solver(x, dim + 1, distance_matrix.row(i_p));
		if (isnan(a[dim]))
			return INFINITY;
		else 
		{
			for (int i_d = 0; i_d < (dim + 1); i_d++)
				x[i_p][i_d] = a[i_d];

			if (a[dim] > 0)
				dim++;
		}

	}

	return dim;
}

double* Solver(double** x, int dim, Eigen::VectorXd distance_vector)
{

	// "distance_vector" has (num_points) elements containing the distance between the (dim) point & the other (dim - 1) points 
	// "x" is (dim, dim-1) giving the coordinates of the (dim) points in the current dimension, (dim-1)
	// So we have (dim-1) equations & (dim-1) unknown coordinates

	double err = 1e-8; // Error Limit

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
		{
			break;
		}

		for (int j = (i + 1); j < dim; j++)
		{
			a[j - 1] = a[j - 1] - x[j][i - 1] * a[i - 1];
		}
	}

	double norm2_a = 0;
	for (int i = 0; i < (dim - 1); i++)
	{
		norm2_a += a[i] * a[i];
	}


	if (abs(norm2_x_new - norm2_a) < err)
	{
		a[dim - 1] = 0.0;
	}
	else if (norm2_x_new > norm2_a)
	{
	a[dim - 1] = sqrt(norm2_x_new - norm2_a);
	}
	else
	{
	a[dim - 1] = NAN;
	}

	return a;
}

double StepFunc_0(double scale, double x)
{
	double out = (double)(x >= scale);
	return out;
}

double StepFunc_1(double scale, double x)
{
	double out;
	if (scale >= 1)
		out = 1;
	else
	{
		if (x <= (1 - scale))
			out = x * scale / (1 - scale);
		else
			out = x * (1 - scale) / scale + (2 * scale - 1) / scale;
	}
	return out;
}

double StepFunc_2(double scale, double x)
{
	double out;
	if (scale <= 0.5)
	{
		if (x <= (1 - 2 * scale) / (1 - scale))
			out = 0;
		else
			out = 1 - (1 - x) * (1 - scale) / scale;
	}
	else
	{
		if (x <= ((1 - scale) / scale))
			out = x * scale / (1 - scale);
		else
			out = 1.0;
	}
	return out;
}

double StepFunc_3(double scale, double x)
{
	double out;
	if (scale <= 0.5)
		out = x * scale / (1 - scale);
	else
		out = 1 - ((1 - x) * (1 - scale) / scale);

	return out;
}

double StepFunc_4(double scale, double x)
{
	double out;
	if ((scale == 1) || (scale == 0))
		out = scale;
	else
	{
		if (scale <= 0.5)
			out = -(1 - x) * (1 - scale) / scale + sqrt((1 - x) * (1 - x) * ((1 - scale) * (1 - scale) / (scale * scale) - 1) + 1);
		else
			out = 1 + x * scale / (1 - scale) - sqrt(x * x * (scale * scale / (1 - scale) / (1 - scale) - 1) + 1);
	}
	return out;
}

double StepFunc_32(double scale, double x)
{
	double out;
	if ((scale > 0.5) && (x > ((1 - scale) / scale)))
		out = 1.0;
	else
		out = x * scale / (1 - scale);

	return out;
}

Multiplicity Compute_Multiplicity(Eigen::VectorXd energy_vector, int num_points, double err)
{
	// * This function assumes that "energy_vector" is sorted *

	Multiplicity M;
	int* initial_ind = new int[num_points];
	int* final_ind = new int[num_points];
	int* nondeg_ind = new int[num_points];
	if (abs(energy_vector(0) - energy_vector(1)) < err) // This is the Initial Index
	{
		M.Num_Degenerate_States++;
		initial_ind[M.Num_Degenerate_States - 1] = 0;
	}
	else // This is a Non-Degenerate Index 
	{
		M.Num_NonDegenerate_States++;
		nondeg_ind[M.Num_NonDegenerate_States - 1] = 0;
	}
	for (int i = 1; i < (num_points - 1); i++)
	{
		if ((abs(energy_vector(i) - energy_vector(i - 1)) > err) && (abs(energy_vector(i) - energy_vector(i + 1)) < err)) // This is the Initial Index 
		{
			M.Num_Degenerate_States++;
			initial_ind[M.Num_Degenerate_States - 1] = i;
		}
		else
		{
			if ((abs(energy_vector(i) - energy_vector(i - 1)) < err) && (abs(energy_vector(i) - energy_vector(i + 1)) > err)) // This is the Final Index 
			{
				final_ind[M.Num_Degenerate_States - 1] = i;
			}
			else if ((abs(energy_vector(i) - energy_vector(i - 1)) > err) && (abs(energy_vector(i) - energy_vector(i + 1)) > err)) // This is a Non-Degenerate Index 
			{
				M.Num_NonDegenerate_States++;
				nondeg_ind[M.Num_NonDegenerate_States - 1] = i;
			}
		}
	}
	if (abs(energy_vector(num_points - 1) - energy_vector(num_points - 2)) < err) // This is the Final Index 
	{
		final_ind[M.Num_Degenerate_States - 1] = num_points - 1;
	}
	else // This is a Non-Degenerate Index 
	{
		M.Num_NonDegenerate_States++;
		nondeg_ind[M.Num_NonDegenerate_States - 1] = num_points - 1;
	}

	M.initial_index = new int[M.Num_Degenerate_States];
	M.final_index = new int[M.Num_Degenerate_States];
	M.degeneracy = new int[M.Num_Degenerate_States];
	for (int i_d = 0; i_d < M.Num_Degenerate_States; i_d++)
	{
		M.initial_index[i_d] = initial_ind[i_d];
		M.final_index[i_d] = final_ind[i_d];
		M.degeneracy[i_d] = final_ind[i_d] - initial_ind[i_d] + 1;
	}

	M.nondegenerate_index = new int[M.Num_NonDegenerate_States];
	for (int i_nd = 0; i_nd < M.Num_NonDegenerate_States; i_nd++)
		M.nondegenerate_index[i_nd] = nondeg_ind[i_nd];


	return M;
}
