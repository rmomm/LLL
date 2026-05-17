#include "LLL.h"
#include <iostream>
#include <vector>
using namespace std;

void printMatrix(const vector<vector<vector<double>>>& mats,
    const string& filename) {

    ofstream out(filename);

    for (size_t k = 0; k < mats.size(); k++) {

        out << "Matrix " << k + 1 << "\n";

        for (const auto& row : mats[k]) {

            for (double x : row) {
                out << x << " ";
            }

            out << "\n";
        }

        out << "\n";
    }

    out.close();
}


int main() {

    auto mats = generateMatrix();

    cout << "Generated matrices: " << mats.size() << endl;

    printMatrix(mats, "matrices.txt");

    cout << "Saved to matrices.txt" << endl;

    return 0;
}