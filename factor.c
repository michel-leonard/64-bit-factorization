// This project is released into the public domain.
// Allowing for unrestricted use, modification, and distribution.

// This file is a C99 implementation of a 64-bit integer factorization algorithm.
// The "factor" algorithm utilizes Pollard's Rho method for factoring large numbers.

typedef unsigned long long int ulong;

typedef struct {
	ulong prime;
	int power;
} row;

ulong mul_mod(ulong a, ulong b, const ulong mod) {
	ulong res = 0, c; // return (a * b) % mod, avoiding overflow errors while doing modular multiplication.
	for (b %= mod; a; a & 1 ? b >= mod - res ? res -= mod : 0, res += b : 0, a >>= 1, (c = b) >= mod - b ? c -= mod : 0, b += c);
	return res % mod;
}

ulong pow_mod(ulong n, ulong exp, const ulong mod) {
	ulong res = 1; // return (n ^ exp) % mod.
	for (n %= mod; exp; exp & 1 ? res = mul_mod(res, n, mod) : 0, n = mul_mod(n, n, mod), exp >>= 1);
	return res;
}

int is_prime(ulong n) {
	// Perform a Miller-Rabin test, it should be a deterministic version.
	static const ulong bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	static const int n_bases = (int) sizeof(*bases);
	for (int i = 0; i < n_bases; ++i)
		if (n % bases[i] == 0)
			return n == bases[i];
	if (n < bases[n_bases - 1] * bases[n_bases - 1])
		return 1 < n;
	// Depending on the size of the number, we don't need to test all the bases.
	int lim = n < 2152302898747 ? n < 25326001 ? n < 2047 ? 1 : n < 1373653 ? 2 : 3 : n < 3215031751 ? 4 :
		5 : n < 341550071728321 ? n < 3474749660383 ? 6 : 7 : n < 3825123056546413051 ? 9 : 12, res = 1, a = 0;
	ulong b, c;
	for (b = n - 1; ~b & 1; b >>= 1, ++a);
	for (int i = 0; i < lim && res; ++i)
		if (c = pow_mod(bases[i], b, n), c != 1) {
			for (int d = a; d-- && (res = c + 1 != n);)
				c = mul_mod(c, c, n);
			res = !res;
		}
	return res;
}

ulong pollard_rho(const ulong n) {
	static const int timeout = 18;
	static ulong r = 88172645463325252; // simple xorshift generator.
	r ^= r << 13, r ^= r >> 7, r ^= r << 17;
	ulong gcd = 1, a, b, c, i = 0, j = 1, x = 1, y = 1 + r % (n - 1);
	for (; gcd == 1; ++i) {
		if (i == j) {
			if (j >> timeout)
				break;
			j <<= 1;
			x = y;
		}
		a = y, b = y;
		for (y = 0; a; a & 1 ? b >= n - y ? y -= n : 0, y += b : 0, a >>= 1, (c = b) >= n - b ? c -= n : 0, b += c);
		y = (1 + y) % n;
		for (a = y > x ? y - x : x - y, b = n; (a %= b) && (b %= a););
		gcd = a | b;
	}
	return gcd;
}

ulong square_root(ulong n) {
	// The Babylonians used an iterative approach for finding square roots.
	ulong a = n != 0, b;
	if (3 < n)
		for (a = n >> 1, b = (a + n / a) >> 1; b < a; a = b, b = (a + n / a) >> 1);
	return a;
}

static inline ulong square_extraction(ulong *n, int *pow) {
	ulong root = 1;
	if (3 < *n)
		while (root = square_root(*n), *n == root * root)
			*n = root, *pow <<= 1;
	return 65522U * 65522U < *n ? 65522 : root + 1;
}

void factor(ulong n, row *rows) {
	if (3 < n) {
		int pow = 1;
		if (~n & 1) {
			// Powers of two are removed.
			*rows = (row) {2, 0};
			do {
				++(*rows).power;
				n >>= 1;
			} while (~n & 1);
			++rows;
		}
		if (8 < n) {
			// The number is odd.
			ulong limit = square_extraction(&n, &pow);
			// Ensure the number has no 16-bit factor.
			for (ulong prime = 3; prime < limit; prime += 2)
				if (n % prime == 0) {
					int p = 0;
					do ++p, n /= prime;
					while (n % prime == 0);
					*rows++ = (row) {prime, p * pow};
					limit = square_extraction(&n, &pow);
				}
			while (n >> 32 && !is_prime(n)) {
				// The number has 2 or 3 prime factors greater than 65536.
				ulong x, y;
				while (x = pollard_rho(n), x == 1 || x == n);
				n /= x;
				if (x >> 32) {
					y = square_root(x);
					if (y * y == x)
						// Pollard's Rho produced a square.
						*rows++ = (row) {y, pow << 1};
					else if (is_prime(x))
						// Pollard's Rho produced a prime.
						*rows++ = (row) {x, pow};
					else {
						// Pollard's Rho produced a composite number.
						while (y = pollard_rho(x), y == 1 || y == x);
						*rows++ = (row) {x / y, pow};
						*rows++ = (row) {y, pow};
					}
				} else if (n % x)
					// Pollard's Rho produced a prime.
					*rows++ = (row) {x, pow};
				else
					// Pollard's Rho produced a prime that divides the number twice.
					n /= x, *rows++ = (row) {x, pow + 1};
			}
		}
		if (n != 1)
			*rows++ = (row) {n, pow};
	} else if (n)
		*rows++ = (row) {n, 1};
	*rows = (row) {0};
}

// Compiled with gcc -Wall -Werror -Wextra -pedantic -O3 -std=c99 factor.c -o fac.exe
// Example of usage :
#include <stdio.h>
int main(void) {
	ulong n = 281496452005891; // a composite number
	row res[16];
	factor(n, res);
	for (int i = 0; res[i].power; ++i)
		if (res[i].power < 2)
			printf("%llu\n", res[i].prime);
		else
			printf("%llu^%d\n", res[i].prime, res[i].power);
	return 0;
}
