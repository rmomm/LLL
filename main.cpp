#include "LLL.h"

#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

void printResults(const vector<vector<vector<double>>>& mats, const string& filename) {
    ofstream out(filename);
    vector<double> deltas ={
        0.5,
        0.75,
        0.90,
        0.95,
        0.99
    };



    for (int i = 0; i < mats.size(); i++) {
        out << "====================================\n";
        out << "Matrix " << i + 1 << "\n\n";

        for (const auto& row : mats[i]) {

            for (double x : row) {

                out << x << " ";
            }

            out << "\n";
        }

        out << "\n";

        for (double delta : deltas) {

            vector<vector<double>> B = mats[i];

            auto start = chrono::high_resolution_clock::now();

            LLLResult res = LLL(B, delta);

            auto end = chrono::high_resolution_clock::now();
            double time_ms = chrono::duration<double, milli> (end - start).count();
            double H = Hadamard(B);
            double b1_norm = vectorNorm(B[0]);


            out << "delta=" << delta << "\n";
            out << "swaps=" << res.swaps << "\n";
            out << "H(B)=" << H << "\n";
            out << "b1_norm=" << b1_norm << "\n";
            out << "time_ms=" << time_ms << "\n\n";
            out << "\n";
        }
    }

    out.close();
}



int main() {
    auto mats = generateMatrix();

    cout << "Generated matrices: " << mats.size() << endl;
    printResults(mats, "results.txt");
    cout << "Saved to results.txt" << endl;

    return 0;
}