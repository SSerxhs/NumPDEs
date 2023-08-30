#pragma once
#include "matrix.h"

valarray<db> restriction(const valarray<db> &a) {
    int n = a.size();
    int N = 1 << (__lg(n) >> 1);
    int M = N >> 1;
    int m = M * M, i, j, x, y;
    assert(N * N == n);
    assert(M * 2 == N);
    valarray<db> b(m);
    for (i = 0; i < M; i++) {
        x = i * N;
        for (j = 0; j < M; j++) {
            y = x + j << 1;
            b[i * M + j] = (a[y] + a[y + 1] + a[y + N] + a[y + N + 1]) / 4;
        }
    }
    return b;
}
valarray<db> interpolation(const valarray<db> &a) {
    int n = a.size();
    int N = 1 << (__lg(n) >> 1);
    assert(N * N == n);
    int M = N << 1;
    int m = M * M, i, j;
    valarray<db> b(m);
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
            b[i * M + j] = a[(i >> 1) * N + (j >> 1)];
    return b;
}
valarray<db> reset_size(valarray<db> a, int N) {
    while (a.size() < N * N)
        a = interpolation(a);
    while (a.size() > N * N)
        a = restriction(a);
    return a;
}
template <typename func_lhs, typename func_rhs>
class multigrid {
    valarray<db> x;
    int times;
    db residual;
    func_lhs get_lhs;
    func_rhs get_rhs;
    void V_cycle(valarray<db> &x, valarray<db> rhs, int v1, int v2) {
        int n = x.size();
        int N = 1 << (__lg(n) >> 1);
        int M = N >> 1;
        int m = M * M, i, j;
        const auto &A = get_lhs(N);
        x = relax(A, x, rhs, v1);
        if (N > 4) {
            valarray rhss = restriction(rhs - multiply(A, x));
            valarray<db> xx(m);
            V_cycle(xx, rhss, v1, v2);
            x += interpolation(xx);
        }
        x = relax(A, x, rhs, v2);
    }
    void FMG(valarray<db> &x, valarray<db> rhs, int v1, int v2) {
        int n = x.size();
        int N = 1 << (__lg(n) >> 1);
        int M = N >> 1;
        int m = M * M, i, j;
        if (N == 4)
            x = 0;
        else {
            valarray rhs2 = restriction(rhs);
            x.resize(m);
            FMG(x, rhs2, v1, v2);
            x = interpolation(x);
        }
        V_cycle(x, rhs, v1, v2);
    }

public:
    multigrid(const func_lhs &_get_lhs, const func_rhs &_get_rhs, const valarray<db> &guess, db eps = 1e-11, int v1 = 4,
        int v2 = 4)
        : get_lhs(_get_lhs), get_rhs(_get_rhs), x(guess) {
        int n = x.size();
        int N = 1 << (__lg(n) >> 1);
        if (N * N != n)
            throw invalid_argument("Incorrect size of guess.");
        int i, j;
        const auto &A = get_lhs(N);
        auto rhs = get_rhs(N);
        rhs -= multiply(A, x);
        for (times = 1; abs(rhs).max() > eps; times++) {
            valarray<db> xx(n);
            FMG(xx, rhs, v1, v2);
            x += xx;
            rhs -= multiply(A, xx);
            // cerr << "current residual = " << abs(rhs).max() << endl;
        }
        residual = abs(rhs).max();
        // cerr << "times = " << times << ", residual = " << residual << endl;
    }
    valarray<db> get_result() const {
        return x;
    }
    int get_times() const {
        return times;
    }
    db get_residual() const {
        return residual;
    }
};

valarray<db> restriction_3d(const valarray<db> &a) {
    int n = a.size();
    int N = 1 << (__lg(n) / 3);
    int M = N >> 1;
    int m = M * M * M, i, j, k, x, y, z;
    // cerr<<"n = "<<n<<endl;
    assert(N * N * N == n);
    assert(M * 2 == N);
    valarray<db> b(m);
    int tp = 0;
    for (i = 0; i < M; i++) {
        x = i * N;
        for (j = 0; j < M; j++) {
            y = (x + j) * N;
            for (k = 0; k < M; k++) {
                z = y + k << 1;
                b[tp++] = (a[z] + a[z + 1] + a[z + N] + a[z + N + 1] + a[z + N * N] + a[z + 1 + N * N] +
                              a[z + N + N * N] + a[z + N + 1 + N * N]) /
                          8;
            }
        }
    }
    return b;
}
valarray<db> interpolation_3d(const valarray<db> &a) {
    int n = a.size();
    int N = 1 << (__lg(n) / 3);
    assert(N * N * N == n);
    int M = N << 1;
    int m = M * M * M, i, j, k, tp = 0;
    valarray<db> b(m);
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++)
            for (k = 0; k < M; k++)
                b[tp++] = a[((i >> 1) * N + (j >> 1)) * N + (k >> 1)];
    return b;
}
valarray<db> reset_size_3d(valarray<db> a, int N) {
    while (a.size() < N * N * N)
        a = interpolation_3d(a);
    while (a.size() > N * N * N)
        a = restriction_3d(a);
    return a;
}
template <typename func_lhs, typename func_rhs>
class multigrid_3d {
    valarray<db> x;
    int times;
    db residual;
    func_lhs get_lhs;
    func_rhs get_rhs;
    void V_cycle(valarray<db> &x, valarray<db> rhs, int v1, int v2) {
        int n = x.size();
        int N = 1 << (__lg(n) / 3);
        int M = N >> 1;
        int m = M * M * M, i, j;
        const auto &A = get_lhs(N);
        x = relax(A, x, rhs, v1);
        if (N > 4) {
            valarray rhss = restriction_3d(rhs - multiply(A, x));
            valarray<db> xx(m);
            V_cycle(xx, rhss, v1, v2);
            x += interpolation_3d(xx);
        }
        x = relax(A, x, rhs, v2);
    }
    void FMG(valarray<db> &x, valarray<db> rhs, int v1, int v2) {
        int n = x.size();
        int N = 1 << (__lg(n) / 3);
        int M = N >> 1;
        int m = M * M * M, i, j;
        if (N == 4)
            x = 0;
        else {
            valarray rhs2 = restriction_3d(rhs);
            x.resize(m);
            FMG(x, rhs2, v1, v2);
            x = interpolation_3d(x);
        }
        V_cycle(x, rhs, v1, v2);
    }

public:
    multigrid_3d(const func_lhs &_get_lhs, const func_rhs &_get_rhs, const valarray<db> &guess, db eps = 1e-11,
        int v1 = 4, int v2 = 4)
        : get_lhs(_get_lhs), get_rhs(_get_rhs), x(guess) {
        int n = x.size();
        int N = 1 << (__lg(n) / 3);
        if (N * N * N != n)
            throw invalid_argument("Incorrect size of guess.");
        int i, j;
        const auto &A = get_lhs(N);
        auto rhs = get_rhs(N);
        rhs -= multiply(A, x);
        for (times = 1; abs(rhs).max() > eps; times++) {
            valarray<db> xx(n);
            FMG(xx, rhs, v1, v2);
            x += xx;
            rhs -= multiply(A, xx);
            // cerr << "current residual = " << abs(rhs).max() << endl;
        }
        residual = abs(rhs).max();
        // cerr << "times = " << times << ", residual = " << residual << endl;
    }
    valarray<db> get_result() const {
        return x;
    }
    int get_times() const {
        return times;
    }
    db get_residual() const {
        return residual;
    }
};
