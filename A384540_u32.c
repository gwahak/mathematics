#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum {
    BLOCK_SIZE = 1 << 16,
    PRIME_LIMIT = 65535,  /* floor(sqrt(UINT32_MAX)) */
    MAX_FACTORS = 9      /* The product of the first ten primes exceeds UINT32_MAX. */
};

typedef struct {
    uint32_t p, e;
} Factor;

typedef struct {
    uint32_t remaining;  /* Zero for inputs with no repeated prime factor. */
    unsigned count;
    Factor factors[MAX_FACTORS];
} Candidate;

static uint32_t primes[PRIME_LIMIT + 1];
static size_t prime_count;

static void build_primes(void)
{
    unsigned char composite[PRIME_LIMIT + 1] = {0};
    for (uint32_t p = 2; p <= PRIME_LIMIT; ++p) {
        if (composite[p])
            continue;
        primes[prime_count++] = p;
        if (p * p <= PRIME_LIMIT) {
            for (uint32_t j = p * p; j <= PRIME_LIMIT; j += p)
                composite[j] = 1;
        }
    }
}

static uint64_t first_multiple(uint32_t lo, uint32_t step)
{
    return ((uint64_t)lo + step - 1) / step * step;
}

static Factor remove_factor(uint32_t *value, uint32_t p)
{
    uint32_t e = 0;
    do {
        *value /= p;
        ++e;
    } while (*value % p == 0);
    return (Factor){p, e};
}

static void factor_block(Candidate *block, uint32_t lo, uint32_t hi)
{
    memset(block, 0, ((size_t)(hi - lo) + 1) * sizeof(*block));

    /* Mark only numbers divisible by a prime square. */
    for (size_t i = 0; i < prime_count; ++i) {
        uint32_t p = primes[i];
        uint32_t square = p * p;
        if (square > hi)
            break;
        for (uint64_t n = first_multiple(lo, square); n <= hi; n += square)
            block[n - lo].remaining = (uint32_t)n;
    }

    /* Sieve in ascending prime order. A factor left over at the end is prime. */
    for (size_t i = 0; i < prime_count; ++i) {
        uint32_t p = primes[i];
        uint32_t square = p * p;
        if (square > hi)
            break;

        uint64_t start = first_multiple(lo, p);
        if (start < square)
            start = square;
        for (uint64_t n = start; n <= hi; n += p) {
            Candidate *candidate = &block[n - lo];
            if (candidate->remaining > 1)
                candidate->factors[candidate->count++] =
                    remove_factor(&candidate->remaining, p);
        }
    }
}

/* Factor value after removing a prime already known to divide it, if supplied. */
static unsigned factorize(uint32_t value, uint32_t known_prime, Factor *factors)
{
    unsigned count = 0;
    if (known_prime)
        factors[count++] = remove_factor(&value, known_prime);
    for (size_t i = 0; i < prime_count; ++i) {
        uint32_t p = primes[i];
        if (p * p > value)
            break;
        if (value % p == 0)
            factors[count++] = remove_factor(&value, p);
    }
    if (value > 1)
        factors[count++] = (Factor){value, 1};
    return count;
}

/* Appending x in ordinary base b scales the prefix by b^digits(x). */
static int append(uint32_t *prefix, uint32_t x, uint32_t b, uint32_t limit)
{
    uint64_t power = b;
    for (uint32_t y = x; y >= b; y /= b)
        power *= b;

    if (power > limit)
        return 0;
    uint64_t value = (uint64_t)*prefix * power + x;
    if (value > limit)
        return 0;
    *prefix = (uint32_t)value;
    return 1;
}

static int matches(uint32_t n, const Candidate *candidate, uint32_t b)
{
    const Factor *f = candidate->factors;
    uint32_t prefix = f[0].p;
    if (f[0].e > 1 && !append(&prefix, f[0].e, b, n))
        return 0;

    for (unsigned i = 1; i < candidate->count; ++i) {
        if (!append(&prefix, f[i].p, b, n))
            return 0;
        if (f[i].e > 1 && !append(&prefix, f[i].e, b, n))
            return 0;
    }
    return prefix == n;
}

/* Each partial product divides n minus the value of the last token. */
static void search_divisors(uint32_t n, const Candidate *candidate,
                            const Factor *factors, unsigned count,
                            uint32_t divisor, uint32_t *best)
{
    if (*best && divisor >= *best)
        return;
    if (count == 0) {
        if (divisor >= 2 && matches(n, candidate, divisor))
            *best = divisor;
        return;
    }

    for (uint32_t e = 0; ; ++e) {
        search_divisors(n, candidate, factors + 1, count - 1, divisor, best);
        if (e == factors->e)
            break;
        divisor *= factors->p;
        if (*best && divisor >= *best)
            break;
    }
}

static uint32_t smallest_base(uint32_t n, Candidate *candidate)
{
    if (candidate->remaining > 1)
        candidate->factors[candidate->count++] = (Factor){candidate->remaining, 1};
    if (candidate->count < 2)
        return 0;

    Factor last = candidate->factors[candidate->count - 1];
    uint32_t last_token = last.e == 1 ? last.p : last.e;

    /* A matching base divides n - last_token. If last.e == 1,
       n - last.p is already known to be divisible by last.p. */
    Factor factors[MAX_FACTORS];
    uint32_t known_prime = last.e == 1 ? last.p : 0;
    unsigned count = factorize(n - last_token, known_prime, factors);
    uint32_t best = 0;
    search_divisors(n, candidate, factors, count, 1, &best);
    return best;
}

static int parse_bound(const char *text, uint32_t *bound)
{
    if (*text < '0' || *text > '9')
        return 0;
    char *end;
    errno = 0;
    uintmax_t value = strtoumax(text, &end, 10);
    if (errno || *end || value == 0 || value > UINT32_MAX)
        return 0;
    *bound = (uint32_t)value;
    return 1;
}

int main(int argc, char **argv)
{
    uint32_t start = 1, limit = UINT32_MAX;
    if (argc > 3 || (argc >= 2 && !parse_bound(argv[1], &start)) ||
        (argc == 3 && !parse_bound(argv[2], &limit)) || start > limit) {
        fprintf(stderr, "Usage: %s [start [limit]], 1 <= start <= limit <= %" PRIu32 "\n",
                argv[0], UINT32_MAX);
        return EXIT_FAILURE;
    }

    Candidate *block = malloc((size_t)BLOCK_SIZE * sizeof(*block));
    if (!block) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    time_t started = time(NULL);
    build_primes();
    uint64_t found = 0;
    const uint64_t progress_step = 10000000;
    uint64_t next_progress = ((uint64_t)start + progress_step - 1)
                           / progress_step * progress_step;

    for (uint64_t lo = start; lo <= limit; lo += BLOCK_SIZE) {
        uint64_t hi = lo + BLOCK_SIZE - 1;
        if (hi > limit)
            hi = limit;

        factor_block(block, (uint32_t)lo, (uint32_t)hi);
        for (uint64_t n = lo; n <= hi; ++n) {
            Candidate *candidate = &block[n - lo];
            if (candidate->remaining) {
                uint32_t b = smallest_base((uint32_t)n, candidate);
                if (b) {
                    printf("%" PRIu64 " : base %" PRIu32 "\n", n, b);
                    fflush(stdout);
                    ++found;
                }
            }
            if (n == next_progress) {
                fprintf(stderr, "Checked %" PRIu64 " (%.1f sec)\n",
                        n, difftime(time(NULL), started));
                next_progress += progress_step;
            }
        }
    }

    printf("Found: %" PRIu64 "\n", found);
    free(block);
    return EXIT_SUCCESS;
}
