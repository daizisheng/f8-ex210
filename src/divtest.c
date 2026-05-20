/* divtest.c -- polynomial divisibility test over F_p, self-contained.
 *
 * Tests whether Q(z)^k divides Q+(z) modulo a prime, by computing the
 * remainder Q+(z) mod Q(z)^k and checking whether it is identically zero.
 *
 * Algorithm is O(deg(Q)^2 * k + deg(Q+) * deg(Q)*k); for the m=5 inputs
 * (deg Q = 8 212, deg Q+ = 39 630, k <= 3) this runs in a few seconds.
 *
 * Compile:   cc -O3 divtest.c -o divtest
 * Usage:     ./divtest <Q.txt> <Q+.txt> <power>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#ifndef MODP
#define MODP 1000000007ULL
#endif

typedef uint64_t u64;

/* Polynomial: coefficient array with logical size n. */

typedef struct { u64 *c; int n; } poly;

static void poly_init(poly *p, int n) {
    p->n = n;
    p->c = calloc((size_t)(n > 0 ? n : 1), sizeof(u64));
}

static void poly_free(poly *p) {
    free(p->c);
    p->c = NULL; p->n = 0;
}

static int poly_degree(const poly *p) {
    for (int i = p->n - 1; i >= 0; i--) if (p->c[i]) return i;
    return -1;
}

/* Modular inverse via Fermat. */

static u64 mod_inv(u64 a, u64 p) {
    u64 result = 1, base = a % p, exp = p - 2;
    while (exp > 0) {
        if (exp & 1) result = result * base % p;
        base = base * base % p;
        exp >>= 1;
    }
    return result;
}

/* Polynomial multiplication mod p (schoolbook). */

static void poly_mul(const poly *a, const poly *b, poly *out) {
    int da = poly_degree(a), db = poly_degree(b);
    if (da < 0 || db < 0) { poly_init(out, 1); return; }
    int n = da + db + 1;
    poly_init(out, n);
    for (int i = 0; i <= da; i++) {
        u64 ai = a->c[i];
        if (!ai) continue;
        for (int j = 0; j <= db; j++) {
            u64 bj = b->c[j];
            if (!bj) continue;
            out->c[i + j] = (out->c[i + j] + ai * bj) % MODP;
        }
    }
}

/* p^k by repeated multiplication. */

static void poly_pow(const poly *a, int k, poly *out) {
    if (k <= 1) {
        poly_init(out, a->n);
        memcpy(out->c, a->c, (size_t)a->n * sizeof(u64));
        return;
    }
    poly tmp;
    poly_init(&tmp, a->n);
    memcpy(tmp.c, a->c, (size_t)a->n * sizeof(u64));
    for (int i = 1; i < k; i++) {
        poly t2;
        poly_mul(&tmp, a, &t2);
        poly_free(&tmp);
        tmp = t2;
    }
    *out = tmp;
}

/* a = b * q + r,  deg r < deg b. */

static void poly_divmod(const poly *a, const poly *b, poly *q, poly *r) {
    int da = poly_degree(a), db = poly_degree(b);
    assert(db >= 0);
    poly_init(r, a->n);
    memcpy(r->c, a->c, (size_t)a->n * sizeof(u64));
    if (da < db) { poly_init(q, 1); return; }
    poly_init(q, da - db + 1);
    u64 lb_inv = mod_inv(b->c[db], MODP);
    for (int i = da - db; i >= 0; i--) {
        u64 coef = r->c[i + db] * lb_inv % MODP;
        q->c[i] = coef;
        if (!coef) continue;
        for (int j = 0; j <= db; j++) {
            u64 sub = coef * b->c[j] % MODP;
            r->c[i + j] = (r->c[i + j] + MODP - sub) % MODP;
        }
    }
}

/* File reader: comments start with #, data lines "index value". */

static int poly_read(const char *path, poly *out) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return -1; }
    int *idx = NULL; u64 *val = NULL;
    int cap = 1024, n = 0, max_idx = -1;
    idx = malloc(cap * sizeof(int));
    val = malloc(cap * sizeof(u64));
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        int i; unsigned long long v;
        if (sscanf(line, "%d %llu", &i, &v) == 2) {
            if (n == cap) {
                cap *= 2;
                idx = realloc(idx, cap * sizeof(int));
                val = realloc(val, cap * sizeof(u64));
            }
            idx[n] = i; val[n] = (u64)v; n++;
            if (i > max_idx) max_idx = i;
        }
    }
    fclose(f);
    poly_init(out, max_idx + 1);
    for (int k = 0; k < n; k++) out->c[idx[k]] = val[k];
    free(idx); free(val);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s <Q.txt> <Q+.txt> <power>\n"
            "Reports whether Q(z)^<power> divides Q+(z) mod %llu.\n",
            argv[0], (unsigned long long)MODP);
        return 1;
    }
    int k = atoi(argv[3]);

    poly Q, Qplus, Q_pow, quot, rem;
    if (poly_read(argv[1], &Q)     < 0) return 1;
    if (poly_read(argv[2], &Qplus) < 0) return 1;

    int dQ     = poly_degree(&Q);
    int dQplus = poly_degree(&Qplus);
    fprintf(stderr, "deg Q = %d, deg Q+ = %d, testing Q^%d divides Q+\n",
            dQ, dQplus, k);

    poly_pow(&Q, k, &Q_pow);
    int dPow = poly_degree(&Q_pow);
    fprintf(stderr, "deg Q^%d = %d\n", k, dPow);

    poly_divmod(&Qplus, &Q_pow, &quot, &rem);
    int dR = poly_degree(&rem);
    if (dR < 0) {
        printf("DIVISIBLE: Q^%d divides Q+\n", k);
    } else {
        printf("NOT DIVISIBLE: Q^%d does not divide Q+ "
               "(remainder degree %d)\n", k, dR);
    }

    poly_free(&Q); poly_free(&Qplus); poly_free(&Q_pow);
    poly_free(&quot); poly_free(&rem);
    return 0;
}
