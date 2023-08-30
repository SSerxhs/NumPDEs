#pragma once
#include "basic.h"
void add(sparse &A, int x, int y, db t) {
    A[x].push_back({y, t});
}
void unique_matrix(sparse &A) {
    for (auto &v : A) {
        vector<pair<int, db>> t;
        swap(v, t);
        int n = t.size(), i, j;
        sort(t.begin(), t.end());
        for (i = j = 0; i < n; i = j) {
            db sum = 0;
            while (j < n && t[i].first == t[j].first)
                sum += t[j++].second;
            if (sum)
                v.push_back({t[i].first, sum});
        }
    }
}
vector<valarray<db>> sparse_to_dense(const sparse &A) {
    int n = A.size(), i;
    vector r(n, valarray<db>(n));
    for (i = 0; i < n; i++)
        for (auto [y, v] : A[i])
            r[i][y] += v;
    return r;
}
sparse dense_to_sparse(const vector<valarray<db>> &A) {
    int n = A.size(), i, j;
    vector<vector<pair<int, db>>> r(n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (A[i][j])
                r[i].push_back({j, A[i][j]});
    return r;
}
valarray<db> relax(const sparse &A, valarray<db> x0, valarray<db> b, int T) {
    int n = x0.size(), i;
    if (n == 0)
        throw invalid_argument("Incorrect length of x0.");
    if (b.size() != n)
        throw invalid_argument("Incorrect length of b.");
    if (A.size() != n)
        throw invalid_argument("Incorrect size of A.");
    valarray<db> D(n);
    for (i = 0; i < n; i++) {
        for (auto [y, v] : A[i]) {
#ifndef NDEBUG
            if (y < 0 || y >= n)
                throw invalid_argument("Incorrect coordinate of A.");
#endif
            if (i == y)
                D[i] += v;
        }
    }
    const static db w = 2.0 / 3, opw = 1 - w;
    auto td = D;
    D = w / D;
    b *= D;
    auto it = A[0].begin();
    while (T--) {
        valarray<db> r(n);
        for (i = 0; i < n; i++)
            for (it = A[i].begin(); it != A[i].end(); ++it)
                r[i] -= it->second * x0[it->first];
        for (i = 0; i < n; i++)
            x0[i] = (r[i] + td[i] * x0[i]) * D[i] + opw * x0[i] + b[i];
    }
    return x0;
}
valarray<db> multiply(const vector<vector<pair<int, db>>> &A, const valarray<db> &x) {
    int n = A.size(), i;
    if (x.size() != n)
        throw invalid_argument("Incorrect length of x.");
    valarray<db> r(n);
    for (i = 0; i < n; i++)
        for (auto [y, v] : A[i]) {
#ifndef NDEBUG0
            if (y < 0 || y >= n)
                throw invalid_argument("Incorrect coordinate of A.");
#endif
            r[i] += x[y] * v;
        }
    return r;
}
