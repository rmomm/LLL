#include "LLL.h"

double dot(const vector<double>& a, const vector<double>& b) {
    double s = 0;
    for (int i = 0; i < a.size(); i++) { 
        s += a[i] * b[i]; 
    }
    return s;
}

double normSq(const vector<double>& v) {
    return dot(v, v);
}

double vectorNorm(const vector<double>& v) {
    return sqrt(normSq(v));
}

void gramSchmidt(const vector<vector<double>>& B, vector<vector<double>>& Bstar, vector<vector<double>>& mu) {
    int m = B.size();      
    int n = B[0].size();

    Bstar.assign(m, vector<double>(n, 0));
    mu.assign(m, vector<double>(m, 0));

    Bstar[0] = B[0];

    for (int i = 1; i < m; i++) {
        Bstar[i] = B[i];

        for (int j = 0; j < i; j++) {
            mu[i][j] = dot(B[i], Bstar[j]) / normSq(Bstar[j]);

            for (int k = 0; k < n; k++) {
                Bstar[i][k] -= mu[i][j] * Bstar[j][k];
            }
        }
    }
}

LLLResult LLL(vector<vector<double>>& B, double delta) {
    int m = B.size();
    int n = B[0].size();
    int swaps = 0;

    vector<vector<double>> Bstar, mu;
    gramSchmidt(B, Bstar, mu);

    int k = 1;

    while (k < m) {
        for (int j = k - 1; j >= 0; j--) {
            if (fabs(mu[k][j]) > 0.5) {
                double q = round(mu[k][j]);

                for (int t = 0; t < n; t++) {
                    B[k][t] -= q * B[j][t];
                }

                gramSchmidt(B, Bstar, mu);
            }
        }


        double left = normSq(Bstar[k]);

        double right = (delta - mu[k][k - 1] * mu[k][k - 1]) * normSq(Bstar[k - 1]);

        if (left >= right) {
            k++;
        }
        else {
            swap(B[k], B[k - 1]);
            swaps++;
            gramSchmidt(B, Bstar, mu);
            k = max(k - 1, 1);
        }
    }

    LLLResult result;
    result.swaps = swaps;
    return result;
}