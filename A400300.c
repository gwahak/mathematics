
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>

#define START 4ULL
#define LIMIT 1073741824ULL

typedef struct {
    uint32_t p, e;
} Factor;

static uint32_t *spf;

/* Smallest prime factor sieve */
void build_sieve(uint32_t limit) {
    spf = calloc((size_t)limit + 1, sizeof(uint32_t));

    if (!spf) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (uint32_t i = 2; i <= limit; i++) {
        if (!spf[i]) {
            spf[i] = i;

            if ((uint64_t)i * i <= limit) {
                for (uint64_t j = (uint64_t)i * i;
                     j <= limit; j += i) {
                    if (!spf[j])
                        spf[j] = i;
                }
            }
        }
    }
}

/* Prime factorization */
int factorize(uint32_t n, Factor *f) {
    int k = 0;

    while (n > 1) {
        uint32_t p = spf[n], e = 0;

        do {
            n /= p;
            e++;
        } while (n > 1 && spf[n] == p);

        f[k++] = (Factor){p, e};
    }

    return k;
}

/*
 * Number of digits of x in bijective base b.
 * Digits range from 1 to b, with no zero.
 */
uint32_t bijective_length(uint32_t x, uint32_t b) {
    uint32_t len = 0;

    do {
        len++;
        x = (x - 1) / b;
    } while (x > 0);

    return len;
}

/*
 * Append the bijective base-b representation of x.
 * The numerical value of the appended block is x.
 */
uint64_t append_bijective(uint64_t v,
                          uint32_t x,
                          uint32_t b) {
    uint32_t len = bijective_length(x, b);
    uint64_t power = 1;

    for (uint32_t i = 0; i < len; i++)
        power *= b;

    return v * power + x;
}

/* Test whether n equals its factorization concatenation */
int matches(uint32_t n, Factor *f, int k, uint32_t b) {
    uint64_t v = 0;

    for (int i = 0; i < k; i++) {
        v = append_bijective(v, f[i].p, b);
        if (v > n) return 0;

        if (f[i].e > 1) {
            v = append_bijective(v, f[i].e, b);
            if (v > n) return 0;
        }
    }

    return v == n;
}

/*
 * Returns the smallest bijective base,
 * or 0 if no base exists.
 */
uint32_t smallest_base(uint32_t n) {
    if (n < 4) return 0;
    if (n == 4) return 1;

    Factor f[32], g[32];
    int k = factorize(n, f);

    /* Exclude primes, but include all composite numbers */
    if (k == 1 && f[0].e == 1)
        return 0;

    uint32_t dr = (f[k-1].e == 1)
                ? f[k-1].p : f[k-1].e;

    if (dr >= n) return 0;

    uint32_t m = n - dr;
    int t = factorize(m, g);

    uint32_t divisors[65536];
    int count = 1;
    divisors[0] = 1;

    for (int i = 0; i < t; i++) {
        int old = count;
        uint32_t power = 1;

        for (uint32_t e = 1; e <= g[i].e; e++) {
            power *= g[i].p;

            for (int j = 0; j < old; j++)
                divisors[count++] = divisors[j] * power;
        }
    }

    uint32_t best = 0;

    for (int i = 0; i < count; i++) {
        uint32_t b = divisors[i];

        if (b >= 2 && (!best || b < best)
            && matches(n, f, k, b))
            best = b;
    }

    return best;
}

int main(void) {
    clock_t start = clock();

    if (LIMIT > UINT32_MAX) {
        fprintf(stderr, "LIMIT exceeds uint32_t\n");
        return 1;
    }

    build_sieve((uint32_t)LIMIT);

    uint64_t found = 0;

    for (uint64_t n = START; n <= LIMIT; n++) {
        uint32_t b = smallest_base((uint32_t)n);

        if (b) {
            printf("%" PRIu64 " : base %" PRIu32 "\n", n, b);
            fflush(stdout);
            found++;
        }

        if (n % 10000000 == 0) {
            fprintf(stderr,
                "Checked %" PRIu64 " (%.1f sec)\n",
                n,
                (double)(clock() - start) / CLOCKS_PER_SEC);
        }
    }

    printf("Found: %" PRIu64 "\n", found);

    free(spf);
    return 0;
}
