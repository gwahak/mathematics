#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN_START (UINT64_C(1) << 32)
#define DEFAULT_LIMIT UINT64_C(13000000000)

enum {
    BLOCK_SIZE = 1 << 16,
    MAX_FACTORS = 15  /* The product of the first 16 primes exceeds UINT64_MAX. */
};

typedef struct {
    uint64_t p;
    unsigned e;
} Factor;

typedef struct {
    uint64_t remaining;  /* Zero for inputs with no repeated prime factor. */
    unsigned count;
    Factor factors[MAX_FACTORS];
} Candidate;

static uint32_t *primes;
static size_t prime_count;

static uint32_t integer_sqrt(uint64_t n)
{
    uint64_t lo = 0, hi = UINT64_C(1) << 32;
    while (hi - lo > 1) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (mid * mid <= n)
            lo = mid;
        else
            hi = mid;
    }
    return (uint32_t)lo;
}

/* Each bit represents one odd integer, with odd n at position n / 2. */
static int is_composite(const uint8_t *sieve, uint64_t n)
{
    return sieve[n / 16] & (1u << ((n / 2) % 8));
}

static int build_primes(uint64_t limit)
{
    uint64_t root = integer_sqrt(limit);
    uint8_t *sieve = calloc((size_t)(root / 16 + 1), sizeof(*sieve));
    if (!sieve)
        return 0;

    for (uint64_t p = 3; p * p <= root; p += 2) {
        if (is_composite(sieve, p))
            continue;
        for (uint64_t n = p * p; n <= root; n += 2 * p)
            sieve[n / 16] |= (uint8_t)(1u << ((n / 2) % 8));
    }

    size_t count = 1;  /* Prime 2. All supported limits have a root above 2. */
    for (uint64_t p = 3; p <= root; p += 2)
        count += (size_t)!is_composite(sieve, p);

    if (count > SIZE_MAX / sizeof(*primes)) {
        free(sieve);
        return 0;
    }
    primes = malloc(count * sizeof(*primes));
    if (!primes) {
        free(sieve);
        return 0;
    }

    primes[prime_count++] = 2;
    for (uint64_t p = 3; p <= root; p += 2) {
        if (!is_composite(sieve, p))
            primes[prime_count++] = (uint32_t)p;
    }
    free(sieve);
    return 1;
}

/* An offset can lie outside the block. Keeping it relative avoids lo + step. */
static uint64_t first_offset(uint64_t lo, uint64_t step)
{
    uint64_t remainder = lo % step;
    return remainder ? step - remainder : 0;
}

static Factor remove_factor(uint64_t *value, uint64_t p)
{
    unsigned e = 0;
    do {
        *value /= p;
        ++e;
    } while (*value % p == 0);
    return (Factor){p, e};
}

static void factor_block(Candidate *block, uint64_t lo, size_t length)
{
    uint64_t hi = lo + (length - 1);
    memset(block, 0, length * sizeof(*block));

    for (size_t i = 0; i < prime_count; ++i) {
        uint64_t p = primes[i];
        uint64_t square = p * p;
        if (square > hi)
            break;
        for (uint64_t j = first_offset(lo, square); j < length;) {
            block[j].remaining = lo + j;
            if (square >= length - j)
                break;
            j += square;
        }
    }

    /* Process primes in ascending order. Any residual factor is prime. */
    for (size_t i = 0; i < prime_count; ++i) {
        uint64_t p = primes[i];
        uint64_t square = p * p;
        if (square > hi)
            break;

        uint64_t offset = square > lo ? square - lo : first_offset(lo, p);
        for (uint64_t j = offset; j < length; j += p) {
            Candidate *candidate = &block[j];
            if (candidate->remaining > 1)
                candidate->factors[candidate->count++] =
                    remove_factor(&candidate->remaining, p);
        }
    }
}

/* The order here is immaterial because these factors only generate divisors. */
static unsigned factorize(uint64_t value, uint64_t known_prime, Factor *factors)
{
    unsigned count = 0;
    if (known_prime)
        factors[count++] = remove_factor(&value, known_prime);
    for (size_t i = 0; i < prime_count; ++i) {
        uint64_t p = primes[i];
        if (p * p > value)
            break;
        if (value % p == 0)
            factors[count++] = remove_factor(&value, p);
    }
    if (value > 1)
        factors[count++] = (Factor){value, 1};
    return count;
}

/* The prefix is positive, so (limit - x) / prefix bounds the base power. */
static int append(uint64_t *prefix, uint64_t x, uint64_t b, uint64_t limit)
{
    if (x > limit)
        return 0;
    uint64_t max_power = (limit - x) / *prefix;
    if (b > max_power)
        return 0;

    uint64_t power = b;
    if (x >= b) {
        uint64_t max_before_multiply = max_power / b;
        for (uint64_t y = x; y >= b; y /= b) {
            if (power > max_before_multiply)
                return 0;
            power *= b;
        }
    }

    *prefix = *prefix * power + x;
    return 1;
}

static int matches(uint64_t n, const Candidate *candidate, uint64_t b)
{
    const Factor *f = candidate->factors;
    uint64_t prefix = f[0].p;
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

/* Every intermediate product divides n - last_token, so it fits in uint64_t. */
static void search_divisors(uint64_t n, const Candidate *candidate,
                            const Factor *factors, unsigned count,
                            uint64_t divisor, uint64_t *best)
{
    if (*best && divisor >= *best)
        return;
    if (count == 0) {
        if (divisor >= 2 && matches(n, candidate, divisor))
            *best = divisor;
        return;
    }

    for (unsigned e = 0; ; ++e) {
        search_divisors(n, candidate, factors + 1, count - 1, divisor, best);
        if (e == factors->e)
            break;
        divisor *= factors->p;
        if (*best && divisor >= *best)
            break;
    }
}

static uint64_t smallest_base(uint64_t n, Candidate *candidate)
{
    if (candidate->remaining > 1)
        candidate->factors[candidate->count++] = (Factor){candidate->remaining, 1};
    if (candidate->count < 2)
        return 0;

    Factor last = candidate->factors[candidate->count - 1];
    uint64_t last_token = last.e == 1 ? last.p : last.e;

    /* A matching base divides n - last_token. If last.e == 1,
       n - last.p is already known to be divisible by last.p. */
    Factor factors[MAX_FACTORS];
    uint64_t known_prime = last.e == 1 ? last.p : 0;
    unsigned count = factorize(n - last_token, known_prime, factors);
    uint64_t best = 0;
    search_divisors(n, candidate, factors, count, 1, &best);
    return best;
}

static int parse_bound(const char *text, uint64_t *bound)
{
    if (*text < '0' || *text > '9')
        return 0;
    char *end;
    errno = 0;
    uintmax_t value = strtoumax(text, &end, 10);
    if (errno || *end || value < MIN_START || value > UINT64_MAX)
        return 0;
    *bound = (uint64_t)value;
    return 1;
}

int main(int argc, char **argv)
{
    uint64_t start = MIN_START, limit = DEFAULT_LIMIT;
    if (argc > 3 || (argc >= 2 && !parse_bound(argv[1], &start)) ||
        (argc == 3 && !parse_bound(argv[2], &limit)) || start > limit) {
        fprintf(stderr, "Usage: %s [start [limit]], %" PRIu64
                        " <= start <= limit <= %" PRIu64 "\n",
                argv[0], MIN_START, UINT64_MAX);
        return EXIT_FAILURE;
    }

    time_t started = time(NULL);
    if (!build_primes(limit)) {
        fprintf(stderr, "Prime table allocation failed\n");
        return EXIT_FAILURE;
    }
    Candidate *block = malloc((size_t)BLOCK_SIZE * sizeof(*block));
    if (!block) {
        fprintf(stderr, "Block allocation failed\n");
        free(primes);
        return EXIT_FAILURE;
    }

    uint64_t found = 0;
    for (uint64_t lo = start; ;) {
        uint64_t remaining = limit - lo;
        size_t length = remaining < BLOCK_SIZE ? (size_t)remaining + 1 : BLOCK_SIZE;
        factor_block(block, lo, length);

        for (size_t j = 0; j < length; ++j) {
            uint64_t n = lo + j;
            Candidate *candidate = &block[j];
            if (candidate->remaining) {
                uint64_t b = smallest_base(n, candidate);
                if (b) {
                    printf("%" PRIu64 " : base %" PRIu64 "\n", n, b);
                    fflush(stdout);
                    ++found;
                }
            }
            if (n % UINT64_C(10000000) == 0) {
                fprintf(stderr, "Checked %" PRIu64 " (%.1f sec)\n",
                        n, difftime(time(NULL), started));
            }
        }
        if (remaining < BLOCK_SIZE)
            break;
        lo += BLOCK_SIZE;
    }

    printf("Found: %" PRIu64 "\n", found);
    free(block);
    free(primes);
    return EXIT_SUCCESS;
}
