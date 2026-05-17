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

struct ExpRes {
    double delta;
    double runtime;
    int swaps;
    double hadamard;
    double b1norm;
};

vector<ExpRes> AllExp(const vector<vector<vector<double>>>& mats) {

    vector<double> deltas =
    { 0.5, 0.75, 0.90, 0.95, 0.99 };

    vector<ExpRes> all;

    for (int i = 0; i < mats.size(); i++) {

        for (double delta : deltas) {

            auto B = mats[i];

            auto start = chrono::high_resolution_clock::now();

            LLLResult res = LLL(B, delta);

            auto end = chrono::high_resolution_clock::now();

            double time_ms = chrono::duration<double, milli>(end - start).count();

            ExpRes e;

            e.delta = delta;
            e.runtime = time_ms;
            e.swaps = res.swaps;
            e.hadamard = Hadamard(B);
            e.b1norm = vectorNorm(B[0]);

            all.push_back(e);
        }
    }

    return all;
}

struct AvgRes {
    double delta;
    double avgTime;
    double avgSwaps;
    double avgHadamard;
    double avgB1;
};

vector<AvgRes> compAvg(const vector<ExpRes>& data) {

    vector<double> deltas = { 0.5, 0.75, 0.90, 0.95, 0.99 };

    vector<AvgRes> result;

    for (double d : deltas) {

        double t = 0, s = 0, h = 0, b = 0;
        int cnt = 0;

        for (auto& x : data) {

            if (fabs(x.delta - d) < 1e-9) {

                t += x.runtime;
                s += x.swaps;
                h += x.hadamard;
                b += x.b1norm;

                cnt++;
            }
        }

        if (cnt == 0) continue;

        AvgRes r;

        r.delta = d;
        r.avgTime = t / cnt;
        r.avgSwaps = s / cnt;
        r.avgHadamard = h / cnt;
        r.avgB1 = b / cnt;

        result.push_back(r);
    }

    return result;
}

void AvgCSV(const vector<AvgRes>& avg, const string& filename) {

    ofstream out(filename);

    out << "delta,avgTime,avgSwaps,avgHadamard,avgB1\n";

    for (auto& r : avg) {

        out << r.delta << ","
            << r.avgTime << ","
            << r.avgSwaps << ","
            << r.avgHadamard << ","
            << r.avgB1 << "\n";
    }

    out.close();
}

int main() {
    auto mats = generateMatrix();

    cout << "Generated matrices: " << mats.size() << endl;
    printResults(mats, "results.txt");
    cout << "Saved to results.txt" << endl;

    auto results = AllExp(mats);
    auto averages = compAvg(results);
    AvgCSV(averages, "averages.csv");

    return 0;
}