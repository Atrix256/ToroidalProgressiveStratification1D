#include "Kritzinger.h"
#include <algorithm>

void KritzingerLDS_AddPoint(std::vector<float>& p, std::vector<float>& sortedp)
{
	// Kritzinger Low Discrepancy Sequence
	// A sequence that outperforms the golden ratio LDS.
	// https://arxiv.org/abs/2406.18132
	//
	// If we have N points in p, we want to find a y value that minimizes the function F(y,p).
	// F(y,p) = (N+1)*y^2 - y - 2*Sum(max(x,y))
	// the Sum is for all x in p.
	// 
	// If we have N points sorted x_1 < ... < x_n, also with x_0 = 0 and x_(n+1) = 1
	// then we have N+1 regions between points.
	// We can solve the function at each region and find whichever is the smallest.
	//
	// To find the minimum of the quadratic, we differentiate, and find where that equals 0.
	//
	// F(y,p) = Ay^2 + By + C
	// F'(y,p) = 2Ay + B
	// 
	// To find the minimum point of F, we find where F' is 0.
	// 2Ay+B = 0 : 2Ay = -B : y = -B/(2A)
	//
	// If y isn't in range of our section, we ignore it.
	// We could clamp it instead of ignoring, but we already have the end points in our point list.
	// If it is in range, we evaluate F at that point and keep the
	// y value that is smallest for all regions.
	// That is our new point

	const size_t N = p.size();

	const float A = float(N + 1);

	float bestY = 0.0f;
	float bestYScore = FLT_MAX;

	for (size_t segmentIndex = 0; segmentIndex <= N; ++segmentIndex)
	{
		// calculate the sum term
		float B = -1.0f;
		float C = 0.0;
		for (size_t xIndex = 0; xIndex < N; ++xIndex)
		{
			if (xIndex < segmentIndex)
				B -= 2.0f;
			else
				C -= 2.0f * sortedp[xIndex];
		}

		// get the edges of our region
		float lastx = (segmentIndex > 0) ? sortedp[segmentIndex - 1] : 0.0f;
		float x = (segmentIndex < N) ? sortedp[segmentIndex] : 1.0f;

		// calculate the minimum y and it's score
		// Verify y is in range
		float y = -B / (2.0f * A);
		if (y <= lastx || y >= x)
			continue;
		float yScore = A * y * y + B * y + C;

		// Keep the best scoring y
		if (yScore < bestYScore)
		{
			bestYScore = yScore;
			bestY = y;
		}
	}

	// Add the point to the list
	p.push_back(bestY);

	// Add the point to the sorted list
	sortedp.push_back(bestY);
	std::sort(sortedp.begin(), sortedp.end());
}