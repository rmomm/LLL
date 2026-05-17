#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

using namespace std;

double dot(const vector<double>& a,const vector<double>& b);
double normSq(const vector<double>& v);
double vectorNorm(const vector<double>& v);
void gramSchmidt( const vector<vector<double>>& B, vector<vector<double>>& Bstar, vector<vector<double>>& mu);


struct LLLResult {
    int swaps;
    double time;
    double hadamard;
    double b1_norm;
};
LLLResult LLL(vector<vector<double>>& B, double delta);


int rankMatrix(vector<vector<double>> A_i);
vector<vector<vector<double>>> generateMatrix();

double determinant(vector<vector<double>> A);
double Hadamard(const vector<vector<double>>& B);
