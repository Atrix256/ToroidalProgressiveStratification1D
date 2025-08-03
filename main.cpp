#define NOMINMAX
#include <Windows.h>

#include <stdio.h>
#include <vector>
#include <string>
#include <random>
#include <direct.h>
#include <unordered_set>
#include <omp.h>
#include <atomic>
#include <numeric>

#include "pcg/pcg_basic.h"

#include "coprime.h"
#include "Kritzinger.h"
#include "ThueMorse.h"

// Random number seed
static const size_t c_randomSeed = 0x1337beef; // 0 for non deterministic
static const bool c_saveMeanAbsError = false;
static const bool c_saveMeanSquaredError = true;
#define MULTI_THREADED() true

using namespace std;

inline float RandomFloat01(pcg32_random_t& rng)
{
	return float(pcg32_random_r(&rng)) / 4294967295.0f;
}

inline pcg32_random_t GetRNG(uint64_t sequence = 0)
{
	pcg32_random_t rng;

	if (c_randomSeed != 0)
	{
		pcg32_srandom_r(&rng, c_randomSeed, sequence);
	}
	else
	{
		std::random_device device;
		std::mt19937 generator(device());
		std::uniform_int_distribution<uint32_t> dist;
		pcg32_srandom_r(&rng, dist(generator), sequence);
	}

	return rng;
}

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

void MakePoints_WhiteNoise(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	for (float& f : points)
		f = RandomFloat01(rng);
}

enum class SequenceOffset
{
	None,
	HalfBucket,
	Random
};

enum class PointOffset
{
	None,
	Stratify
};

enum class Ordering
{
	Sequential,
	ShuffleWhite,
	ShuffleGoldenRatio,
};

class PCGRNG
{
public:
	PCGRNG(pcg32_random_t& rng)
		: m_rng (rng) 
	{
	}

	typedef size_t result_type;
	static size_t min() { return 0; }
	static size_t max() { return 4294967295; }
	size_t operator()()
	{
		return pcg32_random_r(&m_rng);
	}

	pcg32_random_t& m_rng;
};

template <SequenceOffset SequenceOffset, Ordering Ordering, PointOffset PointOffset>
void MakePoints_Regular(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	#if 0
		//std::unordered_set<size_t> indicesSeen;
		//indicesSeen.reserve(numPoints);

		/*
		if (indicesSeen.count(index) > 0)
		{
			printf("Error! already have seen %i\n", (int)index);
		}
		indicesSeen.insert(index);
		*/
	#endif

	pcg32_random_t rng = GetRNG(sequence);

	std::vector<size_t> shuffledPoints;
	if (Ordering == Ordering::ShuffleWhite)
	{
		PCGRNG pcgrng(rng);
		shuffledPoints.resize(numPoints);
		std::iota(shuffledPoints.begin(), shuffledPoints.end(), 0);
		std::shuffle(shuffledPoints.begin(), shuffledPoints.end(), pcgrng);
	}

	size_t grStep = 0;
	if (Ordering == Ordering::ShuffleGoldenRatio)
		grStep = GetIrrationalCoprime(numPoints);

	points.resize(numPoints);

	float sequenceOffset = 0.0f;
	switch (SequenceOffset)
	{
		case SequenceOffset::None: sequenceOffset = 0.0f; break;
		case SequenceOffset::HalfBucket: sequenceOffset = 0.5f / float(numPoints); break;
		case SequenceOffset::Random: sequenceOffset = RandomFloat01(rng); break;
	}

	size_t grIndex = 0;
	for (size_t pointIndex = 0; pointIndex < numPoints; ++pointIndex)
	{
		// Handle point ordering
		size_t locationIndex = 0;
		switch (Ordering)
		{
			case Ordering::Sequential: locationIndex = pointIndex; break;
			case Ordering::ShuffleWhite: locationIndex = shuffledPoints[pointIndex]; break;
			case Ordering::ShuffleGoldenRatio: locationIndex = grIndex; break;
		}
		grIndex = (grIndex + grStep) % numPoints;

		// handle point offsets
		float pointOffset = 0.0f;
		switch (PointOffset)
		{
			case PointOffset::None: pointOffset = 0.0f; break;
			case PointOffset::Stratify: pointOffset = RandomFloat01(rng) / float(numPoints); break;
		}

		// write the point out
		points[pointIndex] = std::fmodf(sequenceOffset + pointOffset + float(locationIndex) / float(numPoints), 1.0f);
	}
}

template <SequenceOffset SequenceOffset, Ordering Ordering, PointOffset PointOffset, int RepeatCount>
void MakePoints_Regular_Repeat(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	points.resize(0);

	for (int iter = 0; iter < RepeatCount; ++iter)
	{
		std::vector<float> newPoints;

		int indexStart = (iter * (int)numPoints / RepeatCount);
		int indexEnd = ((iter + 1) * (int)numPoints / RepeatCount);
		int numIndices = indexEnd - indexStart;

		MakePoints_Regular<SequenceOffset, Ordering, PointOffset>(newPoints, numIndices, sequence * RepeatCount + iter);

		points.insert(points.end(), newPoints.begin(), newPoints.end());
	}
}

static int PermuteNoop(int index, int j, int remainder)
{
	return remainder;
}

template <int Base, typename TPermuteFn>
inline float VanDerCorput(int index, const TPermuteFn& PermuteFn)
{
	float result = 0.0f;
	float f = 1.0f / float(Base);
	float i = float(index);

	for (int x = 0; x < 16; x++)
	{
		if (i <= 0.0f)
			break;

		int remainder = (int)std::fmodf(i, float(Base));

		remainder = PermuteFn(index, x, remainder);

		result += f * std::fmodf(float(remainder), float(Base));
		i = std::floor(i / float(Base));
		f = f / float(Base);
	}

	return result;
}

template <int Base, bool RandomOffset>
void MakePoints_VanDerCorput(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	float offset = RandomOffset ? RandomFloat01(rng) : 0.0f;

	for (size_t i = 0; i < numPoints; ++i)
		points[i] = std::fmodf(offset + VanDerCorput<Base>((int)i, PermuteNoop), 1.0f);
}

// https://perso.liris.cnrs.fr/ostrom/publications/pdf/vo_mcqmc08_article.pdf
template <bool RandomOffset>
void MakePoints_Ostromoukhov60(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	static const int Permutation60[] =
	{
		0,15,30,40,2,48,20,35,8,52,23,43,12,26,55,4,32,45,17,37,
		6,50,28,10,57,21,41,13,33,54,1,25,46,18,38,5,49,29,9,58,
		22,42,14,34,53,3,27,47,16,36,7,51,19,44,31,11,56,24,39,59
	};

	auto Permute =
		[](int index, int j, int remainder)
		{
			return Permutation60[remainder];
		}
	;

	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	float offset = RandomOffset ? RandomFloat01(rng) : 0.0f;

	for (size_t i = 0; i < numPoints; ++i)
		points[i] = std::fmodf(offset + VanDerCorput<60>((int)i, Permute), 1.0f);
}

// https://perso.liris.cnrs.fr/ostrom/publications/pdf/vo_mcqmc08_article.pdf
template <bool RandomOffset>
void MakePoints_Ostromoukhov84(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	static const int Permutation84[] =
	{
		0,22,64,32,50,76,10,38,56,18,72,45,6,28,59,79,41,13,67,25,54,
		2,36,70,16,48,81,30,61,8,43,74,20,52,4,34,66,15,46,77,26,11,62,
		39,82,57,23,69,33,3,51,19,73,42,7,60,29,80,47,14,65,35,1,53,24,
		68,12,40,78,58,27,5,44,71,17,55,37,83,21,49,75,9,31,63
	};

	auto Permute =
		[](int index, int j, int remainder)
		{
			return Permutation84[remainder];
		}
	;

	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	float offset = RandomOffset ? RandomFloat01(rng) : 0.0f;

	for (size_t i = 0; i < numPoints; ++i)
		points[i] = std::fmodf(offset + VanDerCorput<84>((int)i, Permute), 1.0f);
}

template <int Base, bool RandomOffset>
void MakePoints_VDC_ThueMorsePermute(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	float offset = RandomOffset ? RandomFloat01(rng) : 0.0f;

	// For base N, VDC has groups of size N that evenly step between 0 and 1.
	// 
	// Thue Morse base N has has N symbols, and every N symbols contains every symbol [0,N) in a permuted order.
	// TM is "maximally fair" in certain senses because each symbol takes a turn being first.
	//
	// We wil use the N TM symbols for the VDC group to define how the indices are permuted
	//
	auto IndexPermute = [] (int index)
		{
			int ordering = ThueMorse<Base>(index);

			int groupStart = (index / Base) * Base;
			return groupStart + ordering;
		}
	;

	//printf("VDC(%i)\n", Base);

	for (size_t i = 0; i < numPoints; ++i)
	{
		int VDCIndex = IndexPermute((int)i);
		points[i] = std::fmodf(offset + VanDerCorput<Base>(VDCIndex, PermuteNoop), 1.0f);

		//printf("[%i] %i    %0.4f\n", (int)i, VDCIndex, points[i]);
	}
	//printf("\n\n");
}

// Make the points once and then re-use them, since they are a bit costly to create
std::vector<float> g_KritzingerCached;
void InitKritzinger(size_t numPoints)
{
	if (g_KritzingerCached.size() == numPoints)
		return;

	printf("Making Kritzinger:\n");
	g_KritzingerCached.clear();
	g_KritzingerCached.reserve(numPoints);
	std::vector<float> sortedPoints;
	int lastPercent = -1;
	for (size_t i = 0; i < numPoints; ++i)
	{
		int percent = int(100.0f * float(i) / float(numPoints));
		if (lastPercent != percent)
		{
			printf("\r%i%%", percent);
			lastPercent = percent;
		}
		KritzingerLDS_AddPoint(g_KritzingerCached, sortedPoints);
	}
	printf("\r100%%\n\n");
}

template <bool RandomOffset>
void MakePoints_Kritzinger(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	points = g_KritzingerCached;

	// Apply a random offset if we should
	if (RandomOffset)
	{
		pcg32_random_t rng = GetRNG(sequence);
		float offset = RandomFloat01(rng);

		for (float& f : points)
			f = std::fmodf(f + offset, 1.0f);
	}
}

template <bool RandomOffset>
void MakePoints_GoldenRatio(std::vector<float>& points, size_t numPoints, uint64_t sequence)
{
	pcg32_random_t rng = GetRNG(sequence);

	points.resize(numPoints);

	float offset = RandomOffset ? RandomFloat01(rng) : 0.0f;

	float value = offset;

	for (float& f : points)
	{
		f = value;
		value = std::fmodf(value + c_goldenRatio, 1.0f);
	}
}

static const float c_actualValue_Triangle = 0.5f;
float Function_Triangle(float x)
{
	return x;
}

static const float c_actualValue_Step = 0.4f;
float Function_Step(float x)
{
	return (x < 0.4) ? 1.0f : 0.0f;
}

static const float c_actualValue_Sine = std::sinf(2.0f) * std::sinf(2.0f) / 2.0f;
float Function_Sine(float x)
{
	return std::sin(4.0f * x);
}

static const float c_actualValue_Gauss = 0.313309f;
float Function_Gauss(float x)
{
	x -= 0.5f;
	return std::exp(-32.0f * x * x);
}

struct TestResults
{
	string name;
	vector<float> meanAbsError;
	vector<float> meanSquareError;
};

template <typename TMakePointsFn, typename TFunctionFn>
TestResults DoTest(const char* name, TMakePointsFn& MakePoints, TFunctionFn& Function, float actualFunctionValue, size_t numPoints, size_t numTests)
{
	printf("%s\n", name);

	TestResults ret;
	ret.name = name;
	ret.meanAbsError.resize(numPoints);
	ret.meanSquareError.resize(numPoints);

	// Do each test in parallel
	std::vector<std::vector<float>> testAbsError(numTests);
	std::vector<std::vector<float>> testSquareError(numTests);
	{
		std::atomic<int> testsDone = 0;
		#if MULTI_THREADED()
		#pragma omp parallel for
		#endif
		for (int testIndex = 0; testIndex < (int)numTests; ++testIndex)
		{
			// report progress
			int lastTestsDone = testsDone.fetch_add(1);
			int lastPercent = int(100.0f * float(lastTestsDone) / float(numTests));
			int percent = int(100.0f * float(lastTestsDone + 1) / float(numTests));
			if (lastPercent != percent)
				printf("\r%i%%", percent);

			// Make the points
			vector<float> points(numPoints);
			MakePoints(points, numPoints, testIndex);

			// run the function and get error
			std::vector<float>& absError = testAbsError[testIndex];
			absError.resize(numPoints);
			std::vector<float>& squareError = testSquareError[testIndex];
			squareError.resize(numPoints);

			float integration = 0.0f;
			for (size_t pointIndex = 0; pointIndex < numPoints; ++pointIndex)
			{
				float y = Function(points[pointIndex]);
				integration = Lerp(integration, y, 1.0f / float(pointIndex + 1));
				float error = integration - actualFunctionValue;
				absError[pointIndex] = std::abs(error);
				squareError[pointIndex] = error * error;
			}
		}

		printf("\r100%%\n");
	}

	// average across tests
	{
		std::atomic<int> testsDone = 0;
		#if MULTI_THREADED()
		#pragma omp parallel for
		#endif
		for (int pointIndex = 0; pointIndex < (int)numPoints; ++pointIndex)
		{
			// report progress
			int lastTestsDone = testsDone.fetch_add(1);
			int lastPercent = int(100.0f * float(lastTestsDone) / float(numPoints));
			int percent = int(100.0f * float(lastTestsDone + 1) / float(numPoints));
			if (lastPercent != percent)
				printf("\r%i%%", percent);

			// average across tests
			float meanAbsError = 0.0f;
			float meanSquareError = 0.0f;
			for (int testIndex = 0; testIndex < (int)numTests; ++testIndex)
			{
				meanAbsError = Lerp(meanAbsError, testAbsError[testIndex][pointIndex], 1.0f / float(testIndex + 1));
				meanSquareError = Lerp(meanSquareError, testSquareError[testIndex][pointIndex], 1.0f / float(testIndex + 1));
			}
			ret.meanAbsError[pointIndex] = meanAbsError;
			ret.meanSquareError[pointIndex] = meanSquareError;
		}

		printf("\r100%%\n");
	}

	return ret;
}

enum class DataSource
{
	MeanAbsError,
	MeanSquaredError
};

void SaveCSV(const char* fileName, const std::vector<TestResults>& testResults, bool writeOORootN, DataSource dataSource)
{
	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	// write the header
	{
		fprintf(file, "\"Index\"");

		if (writeOORootN)
			fprintf(file, ",\"OneOverRootN\"");

		for (const TestResults& results : testResults)
			fprintf(file, ",\"%s\"", results.name.c_str() );
		fprintf(file, "\r\n");
	}

	// Write each line of data
	for (size_t row = 0; row < testResults[0].meanAbsError.size(); ++row)
	{
		fprintf(file, "\"%zi\"", row + 1);
		if (writeOORootN)
			fprintf(file, ",\"%f\"", 1.0f / std::sqrt(float(row + 1)));
		for (const TestResults& results : testResults)
			fprintf(file, ",\"%f\"", (dataSource == DataSource::MeanAbsError) ? results.meanAbsError[row] : results.meanSquareError[row]);
		fprintf(file, "\r\n");
	}

	fclose(file);
}

void RunCommandLine(const char* commandLine)
{
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;

	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));

	if (!CreateProcessA(
		NULL,                 // No module name (use command line)
		(LPSTR)commandLine,   // Command line
		NULL,                 // Process handle not inheritable
		NULL,                 // Thread handle not inheritable
		FALSE,                // Set handle inheritance to FALSE
		CREATE_NO_WINDOW,     // Creation flags (e.g., no window)
		NULL,                 // Use parent's environment block
		NULL,                 // Use parent's starting directory
		&si,                  // Pointer to STARTUPINFO structure
		&pi                   // Pointer to PROCESS_INFORMATION structure
	)) {
		return;
	}

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// Toroidally progressive stratification tests
	{
		static const size_t numPoints = 100;
		static const size_t numTests = 10000;

		struct Test
		{
			float (*function)(float) = nullptr;
			const char* functionName = nullptr;
			const char* equation = nullptr;
			float actualValue = 0.0f;
		};

		static const Test tests[] =
		{
			{ Function_Triangle, "Triangle", "y = x\nx in [0,1)", c_actualValue_Triangle},
			{ Function_Step, "Step", "y = (x < 0.4) ? 1.0 : 0.0\nx in [0,1)", c_actualValue_Step},
			{ Function_Sine, "Sine", "y = sin(4x)\nx in [0,1)", c_actualValue_Sine},
			{ Function_Gauss, "Gauss", "y = e^(-32(x-0.5)^2)\nx in [0,1)", c_actualValue_Gauss},
		};

		for (const Test& test : tests)
		{
			printf("[%s]\nDoing %i tests for %i points:\n\n", test.functionName, (int)numTests, (int)numPoints);
			std::vector<TestResults> results;
			results.push_back(DoTest("WhiteNoise", MakePoints_WhiteNoise, test.function, test.actualValue, numPoints, numTests));

			results.push_back(DoTest("GoldenRatio", MakePoints_GoldenRatio<true>, test.function, test.actualValue, numPoints, numTests));
			
			//InitKritzinger(numPoints);
			//results.push_back(DoTest("Kritzinger", MakePoints_Kritzinger<true>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("Ostromoukhov60", MakePoints_Ostromoukhov60<true>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("Ostromoukhov84", MakePoints_Ostromoukhov84<true>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("RegularZero", MakePoints_Regular<SequenceOffset::None, Ordering::Sequential, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("Regular", MakePoints_Regular<SequenceOffset::HalfBucket, Ordering::Sequential, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("RegularRandOffset", MakePoints_Regular<SequenceOffset::Random, Ordering::Sequential, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("Stratification", MakePoints_Regular<SequenceOffset::None, Ordering::Sequential, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));
			results.push_back(DoTest("Stratified", MakePoints_Regular<SequenceOffset::Random, Ordering::Sequential, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));//results.push_back(DoTest("StratificationRandOffset", MakePoints_Regular<SequenceOffset::Random, Ordering::Sequential, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("RegularSequential", MakePoints_Regular<SequenceOffset::None, Ordering::Sequential, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("RegularWhite", MakePoints_Regular<SequenceOffset::None, Ordering::ShuffleWhite, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("RegularGR", MakePoints_Regular<SequenceOffset::None, Ordering::ShuffleGoldenRatio, PointOffset::None>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("StratifySequential", MakePoints_Regular<SequenceOffset::Random, Ordering::Sequential, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("StratifyWhite", MakePoints_Regular<SequenceOffset::Random, Ordering::ShuffleWhite, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));
			results.push_back(DoTest("StratifiedGR", MakePoints_Regular<SequenceOffset::Random, Ordering::ShuffleGoldenRatio, PointOffset::Stratify>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("StratifySequentialRepeat2", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::Sequential, PointOffset::Stratify, 2>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("StratifyWhiteRepeat2", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::ShuffleWhite, PointOffset::Stratify, 2>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("StratifiedGRRepeat2", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::ShuffleGoldenRatio, PointOffset::Stratify, 2>, test.function, test.actualValue, numPoints, numTests));


			//results.push_back(DoTest("StratifySequentialRepeat10", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::Sequential, PointOffset::Stratify, 10>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("StratifyWhiteRepeat10", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::ShuffleWhite, PointOffset::Stratify, 10>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("StratifiedGRRepeat10", MakePoints_Regular_Repeat<SequenceOffset::Random, Ordering::ShuffleGoldenRatio, PointOffset::Stratify, 10>, test.function, test.actualValue, numPoints, numTests));
			

			//results.push_back(DoTest("VDC2", MakePoints_VanDerCorput<2, true>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("VDC3", MakePoints_VanDerCorput<3, true>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("VDC7", MakePoints_VanDerCorput<7, true>, test.function, test.actualValue, numPoints, numTests));

			//results.push_back(DoTest("TMVDC2", MakePoints_VDC_ThueMorsePermute<2, true>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("TMVDC3", MakePoints_VDC_ThueMorsePermute<3, true>, test.function, test.actualValue, numPoints, numTests));
			//results.push_back(DoTest("TMVDC7", MakePoints_VDC_ThueMorsePermute<7, true>, test.function, test.actualValue, numPoints, numTests));

			printf("\nSaving csvs and making graphs\n");

			// Mean Abs Error
			if (c_saveMeanAbsError)
			{
				char fileName[1024];
				sprintf_s(fileName, "out/%i_%s.meanAbsError.csv", (int)numPoints, test.functionName);
				SaveCSV(fileName, results, false, DataSource::MeanAbsError);

				printf("Making Graph\n\n");
				char buffer[1024];
				sprintf_s(buffer, "python csvlogloggraph.py %s \"%s: Mean Abs Error Over %i Tests \n%s\"", fileName, test.functionName, (int)numTests, test.equation);
				RunCommandLine(buffer);
			}

			// Mean Squared Error
			if (c_saveMeanSquaredError)
			{
				char fileName[1024];
				sprintf_s(fileName, "out/%i_%s.meanSquaredError.csv", (int)numPoints, test.functionName);
				SaveCSV(fileName, results, false, DataSource::MeanSquaredError);

				printf("Making Graph\n\n");
				char buffer[1024];
				sprintf_s(buffer, "python csvlogloggraph.py %s \"%s: Mean Squared Error Over %i Tests \n%s\"", fileName, test.functionName, (int)numTests, test.equation);
				RunCommandLine(buffer);
			}
		}
	}

	return 0;
}
