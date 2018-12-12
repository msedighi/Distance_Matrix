#pragma once
#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H
#endif // !1

#include "Eigen/Dense"
#include "Math.h"

typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;
typedef Eigen::Array<bool, 1, Eigen::Dynamic> ArrayXb;

class Dendo
{
private:
	int Num_Points;
public:
	bool** Structure;
	double* Scale;

	Dendo(int num_points);
	~Dendo();
};
struct Distance_Index
{
	double Distance;
	unsigned int Index1;
	unsigned int Index2;
};
struct Distance_Struct
{
	Eigen::MatrixXd Operator;
	Distance_Index * Vector;

	long Number_Points;
	long Number_Pairs;

	Distance_Struct(long N)
	{
		Number_Points = N;
		Number_Pairs = N * (N - 1) / 2;

		Operator = Eigen::MatrixXd::Constant(Number_Points, Number_Points, 0);
		Vector = new Distance_Index[Number_Pairs];
	}
};

//
double Euclidean_Distance(double* v1, double* v2, int dim);

double Euclidean_Distance_Vector(Eigen::VectorXd v1, Eigen::VectorXd v2);

template<int p>
double lpDistance(Eigen::VectorXd v1, Eigen::VectorXd v2);

int compareDistance(const void * a, const void * b);

template<class MyType>
int compareMyType(const void * a, const void * b);
//

double*** Difference_Matrix(double** vectors, int num_points, int dim);

double* Vector_Product(double** x, double** y, int num_points, int dim);

Distance_Struct Distance_Matrix(double(*Distance_func)(double*, double*, int), double** Positions, long Number_Points, int dim);

int Embedding_Dimension(Eigen::MatrixXd distance_matrix, int num_points, double** x);

double* Solver(double** x, int dim, Eigen::VectorXd distance_vector);

double StepFunc_0(double scale, double x); // Original Step Function
double StepFunc_1(double scale, double x); // 2 Lines
double StepFunc_2(double scale, double x); // Constant + Line 
double StepFunc_3(double scale, double x); // Linear
double StepFunc_4(double scale, double x); // Heyperbolic 

double StepFunc_32(double scale, double x);

struct Multiplicity
{
	int Num_Degenerate_States = 0;
	int* initial_index;
	int* final_index;
	int* degeneracy;

	int Num_NonDegenerate_States = 0;
	int* nondegenerate_index;
};
Multiplicity Compute_Multiplicity(Eigen::VectorXd energy_vector, int num_points, double err);