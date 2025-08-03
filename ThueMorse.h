#pragma once

// returns the Thue Morse value at the given index, for the given base.
// In base 2, the Thue Morse sequence starts: 01101001...
// In the above, sequence[3] is 0, so the function call ThueMorse<2>(3) would return 0.
// It does this by converting the index to the given base, summing up the digits, and returning that sum, modulo the base.
template <int Base>
int ThueMorse(int index)
{
	int sum = 0;
	while (index)
	{
		sum += index % Base;
		index /= Base;
	}
	return sum % Base;
}
