#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef START
#define START 4ULL
#endif

#ifndef LIMIT
#define LIMIT UINT32_MAX
#endif

typedef struct {
    uint32_t p, e;
} Factor;

static uint16_t *spf;

/*
 * Only numbers congruent to 1 or 5 modulo 6 need a table entry.
 * Their n / 3 indices are consecutive. Zero denotes an unmarked number.
 * Every composite uint32_t has a smallest prime factor <= 65535.
 */
static void build_sieve(uint32_t limit) {
    const size_t entries = (size_t)(limit / 3) + 1;
    spf = calloc(entries, sizeof(*spf));

    if (!spf) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (uint32_t p = 5, prime_step = 2;
         (uint64_t)p * p <= limit; p += prime_step, prime_step = 6 - prime_step) {
        if (spf[p / 3])
            continue;

        /*
         * Start at p*p and visit p times each subsequent cofactor
         * congruent to 1 or 5 modulo 6.
         */
        const uint32_t period = 6 * p;
        uint32_t step = p * (prime_step == 2 ? 2 : 4);

        for (uint64_t j = (uint64_t)p * p; j <= limit;
             j += step, step = period - step) {
            uint16_t *entry = &spf[j / 3];
            if (!*entry)
                *entry = (uint16_t)p;
        }
    }
}

static uint32_t least_factor(uint32_t n) {
    if (n % 2 == 0)
        return 2;
    if (n % 3 == 0)
        return 3;

    const uint32_t p = spf[n / 3];
    return p ? p : n;
}

static int factorize(uint32_t n, Factor *f) {
    int k = 0;

    while (n > 1) {
        const uint32_t p = least_factor(n);
        uint32_t e = 0;

        do {
            n /= p;
            e++;
        } while (n > 1 && n % p == 0);

        f[k++] = (Factor){p, e};
    }

    return k;
}

/*
 * Append the bijective base-b representation of x.
 * Return n + 1 if the new value necessarily exceeds n.
 */
static uint64_t append_bijective(uint64_t v, uint32_t x,
                                 uint32_t b, uint32_t n) {
    if (v == 0)
        return x;

    uint64_t power = 1;
    uint32_t remaining = x;

    do {
        power *= b;
        if (power > n)
            return (uint64_t)n + 1;

        remaining = (remaining - 1) / b;
    } while (remaining > 0);

    /* Both v and power are at most UINT32_MAX, and x is uint32_t. */
    return v * power + x;
}

static int matches(uint32_t n, const Factor *f, int k, uint32_t b) {
    uint64_t v = 0;

    for (int i = 0; i < k; i++) {
        v = append_bijective(v, f[i].p, b, n);
        if (v > n)
            return 0;

        if (f[i].e > 1) {
            v = append_bijective(v, f[i].e, b, n);
            if (v > n)
                return 0;
        }
    }

    return v == n;
}

static uint32_t smallest_base(uint32_t n) {
    if (n < 4)
        return 0;
    if (n == 4)
        return 1;
    if (least_factor(n) == n)
        return 0;

    Factor f[32], g[32];
    const int k = factorize(n, f);
    const Factor last = f[k - 1];
    const uint32_t dr = last.e == 1 ? last.p : last.e;

    if (dr >= n)
        return 0;

    /* The final block is dr modulo b, so every candidate b divides n - dr. */
    const int t = factorize(n - dr, g);
    uint32_t divisors[65536];
    size_t count = 1;
    divisors[0] = 1;

    for (int i = 0; i < t; i++) {
        const size_t old = count;
        uint32_t power = 1;

        for (uint32_t e = 1; e <= g[i].e; e++) {
            power *= g[i].p;

            for (size_t j = 0; j < old; j++)
                divisors[count++] = divisors[j] * power;
        }
    }

    uint32_t best = 0;

    for (size_t i = 1; i < count; i++) {
        const uint32_t b = divisors[i];

        if ((!best || b < best) && matches(n, f, k, b))
            best = b;
    }

    return best;
}

int main(void) {
    const clock_t start = clock();

    if ((uint64_t)START > (uint64_t)LIMIT ||
        (uint64_t)LIMIT > UINT32_MAX) {
        fprintf(stderr, "Require 0 <= START <= LIMIT <= UINT32_MAX\n");
        return EXIT_FAILURE;
    }

    build_sieve((uint32_t)LIMIT);

    uint64_t found = 0;

    for (uint64_t n = START; n <= (uint64_t)LIMIT; n++) {
        const uint32_t b = smallest_base((uint32_t)n);

        if (b) {
            printf("%" PRIu64 " : base %" PRIu32 "\n", n, b);
            fflush(stdout);
            found++;
        }

        if (n % 10000000 == 0) {
            fprintf(stderr, "Checked %" PRIu64 " (%.1f sec)\n",
                    n, (double)(clock() - start) / CLOCKS_PER_SEC);
        }
    }

    printf("Found: %" PRIu64 "\n", found);
    free(spf);
    return EXIT_SUCCESS;
}
