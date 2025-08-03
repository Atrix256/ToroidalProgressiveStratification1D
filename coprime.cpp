#include "coprime.h"

#include <algorithm>

// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
static inline size_t CalculateGCD(size_t smaller, size_t larger)
{
	// make sure A <= B before starting
	if (larger < smaller)
		std::swap(smaller, larger);

	// loop
	while (1)
	{
		// if the remainder of larger / smaller is 0, they are the same
		// so return smaller as the GCD
		size_t remainder = larger % smaller;
		if (remainder == 0)
			return smaller;

		// otherwise, the new larger number is the old smaller number, and
		// the new smaller number is the remainder
		larger = smaller;
		smaller = remainder;
	}
}

static inline bool IsCoprime(size_t A, size_t B)
{
	return CalculateGCD(A, B) == 1;
}

size_t GetIrrationalCoprime(size_t n, float irrational)
{
	size_t target = size_t(std::fmod(irrational, 1.0f) * float(n) + 0.5f);

	size_t offset = 0;
	while (1)
	{
		if (offset < target)
		{
			size_t coprime = target - offset;
			if (IsCoprime(coprime, n))
				return coprime;
		}

		size_t coprime = target + offset + 1;
		if (coprime < n && IsCoprime(coprime, n))
			return coprime;

		offset++;
	}
}