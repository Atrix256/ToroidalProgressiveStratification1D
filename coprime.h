#pragma once

#include <cmath>

#define c_goldenRatio ((1.0f + std::sqrt(5.0f)) / 2.0f)

// returns a number in the ring Zn [0,n) that is coprime to n, closest to the irrational given if you were to map [0, n) to [0, 1).
size_t GetIrrationalCoprime(size_t n, float irrational = c_goldenRatio);
