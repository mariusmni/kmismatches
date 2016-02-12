#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

#include "kutil.h"
#include "algnaive.h"
#include "algabrahamson.h"
#include "algnk.h"
#include "algnsqrtk.h"
#include "algrandomized.h"

int K[] = { 1, 5, 10, 50};

void readInput(char *inFile, int *&t, int &n, int *&p, int &m, int &k) {
	FILE *f = fopen(inFile, "r");
	if (f == NULL) {
		cout << "Cannot open file " << inFile << endl;
		exit(0);
	}
	fscanf(f, "%d %d %d\n", &m, &n, &k);
	p = new int[m + 1];
	for (int i = 0; i < m; ++i)
		fscanf(f, "%d", p + i);
	t = new int[n + 1];
	for (int i = 0; i < n; ++i)
		fscanf(f, "%d", t + i);
	fclose(f);
}

void output(FILE *f, int *a, int n) {
	for (int i = 0; i < n; ++i)
		fprintf(f, "%d ", a[i]);
	fprintf(f, "\n");
}

double frand(void) {
	double value;
	value = ((double) rand() / (RAND_MAX));
	return value;
}

void convertTo01(int *m, int n, int k) {
	for (int i = 0; i < n; ++i)
		m[i] = m[i] <= k ? 1 : 0;
}

void trim(int *m, int n, int k) {
	for (int i = 0; i < n; ++i)
		if (m[i] > k)
			m[i] = k;
}

void test(int *t, int n, int *p, int m, int k) {
	cout << n << " \t " << m << " \t " << k;
	resetTime();
	int *m1 = naiveCountMismatches(t, n, p, m);
	cout << " \t " << getTime();
//		printTime("Naive");

//	if (false)
	{// Abrahamson
		resetTime();
		int *m2 = abrahamsonCountMismatches(t, n, p, m);
		cout << " \t " << getTime();
//			printTime("Abrahamson");
		if (memcmp(m1, m2, n * sizeof(int)) != 0) {
			cout << "Error in abrahamson" << endl;
			int i = 0;
			while (m2[i] == m1[i])
				++i;
			print("Text      ", t, n);
			print("Pattern   ", p, m);
			print("Naive     ", m1, n);
			print("Abrah     ", m2, n);
			cout << "First difference at " << i << endl;
			print("Text diff ", t + i, n - i);
			print("Pattern   ", p, m);
			exit(1);
		}
		delete[] m2;
	}

	trim(m1, n, k + 1);

//	if (false)
	{ // O(nk) with m^2 preprocessing
		resetTime();
		int *flag = new int[n];
		for (int i = 0; i < n; ++i)
			flag[i] = 1;
		m2preproc m2pp(p, m);
		int *m4 = nkCountMismatches(t, n, p, m, k, flag, m2pp);
		cout << " \t " << getTime();
//			printTime("NK m2preproc");
		if (memcmp(m4, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nk with m2preproc and naive" << endl;
			int i = 0;
			while (m4[i] == m1[i])
				++i;
			cout << "First difference at " << i << " naive said " << m1[i]
					<< " other " << m4[i] << endl;
			print("Text diff ", t + i, min(m, n - i));
			print("Pattern   ", p, m);
			//			exit(1);
		}
		delete[] m4;
		delete[] flag;
	}

//	if (false)
	{ // O(nk) with suffix array preprocessing
		resetTime();
		int *flag = new int[n];
		for (int i = 0; i < n; ++i)
			flag[i] = 1;
		sapreproc sapp(p, m);
		int *m5 = nkCountMismatches(t, n, p, m, k, flag, sapp);
		cout << " \t " << getTime();
//			printTime("NK sapreproc");
		if (memcmp(m5, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nk with sapreproc and naive" << endl;
			exit(1);
		}
		delete[] m5;
		delete[] flag;
	}

	if (false)
	{ // O(nk) with suffix tree preprocessing
		resetTime();
		int *flag = new int[n];
		for (int i = 0; i < n; ++i)
			flag[i] = 1;
		stpreproc stpp(p, m);
		int *m8 = nkCountMismatches(t, n, p, m, k, flag, stpp);
		cout << " \t " << getTime();
		//	printTime("NK sapreproc");
		if (memcmp(m8, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nk with stpreproc and naive" << endl;
			int i = 0;
			while (m8[i] == m1[i])
				++i;
			cout << "First difference at " << i << " naive said " << m1[i]
					<< " other " << m8[i] << endl;
			print("Text ", t, n);//i, min(m, n - i));
			print("Pattern   ", p, m);
			exit(1);
		}
		delete[] m8;
		delete[] flag;
	}

	{ // O(nsqrtk) with O(m^2) preprocessing
		resetTime();
		int *m6 = nsqrtKMismatches(t, n, p, m, k, 0);
		cout << " \t " << getTime();
		trim(m6, n, k + 1);
//			printTime("NsqrtK m2preproc");
		if (memcmp(m6, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nsqrtk and naive" << endl;
			int i = 0;
			while (m6[i] == m1[i])
				++i;
			//			print("Text      ", t, n);
			//			print("Pattern   ", p, m);
			//			print("Naive     ", m1, n);
			//			print("NsqrtK    ", m6, n);
			cout << "First difference at " << i << "naive said" << m1[i]
					<< " other " << m6[i] << endl;
			print("Text diff ", t + i, n - i);
			print("Pattern   ", p, m);
			exit(1);
		}
		delete[] m6;
	}

	{ // O(nsqrtk) with suffix array preprocessing
		resetTime();
		int *m7 = nsqrtKMismatches(t, n, p, m, k, 1);
		cout << " \t " << getTime();
		trim(m7, n, k + 1);
		//	printTime("NsqrtK sapreproc");
		if (memcmp(m7, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nsqrtk sa and naive" << endl;
			int i = 0;
			while (m7[i] == m1[i])
				++i;
			print("Text      ", t, n);
			print("Pattern   ", p, m);
			print("Naive     ", m1, n);
			print("NsqrtK    ", m7, n);
			cout << "First difference at " << i << "naive said" << m1[i]
					<< " other " << m7[i] << endl;
			print("Text diff ", t + i, n - i);
			print("Pattern   ", p, m);
			exit(1);
		}
		delete[] m7;
	}

	if (false)
	{ // O(nsqrtk) with suffix tree preprocessing
		resetTime();
		int *m9 = nsqrtKMismatches(t, n, p, m, k, 2);
		cout << " \t " << getTime();
		trim(m9, n, k + 1);
			printTime("NsqrtK sapreproc");
		if (memcmp(m9, m1, n * sizeof(int)) != 0) {
			cout << "Difference between nsqrtk st and naive" << endl;
			int i = 0;
			while (m9[i] == m1[i])
				++i;
			print("Text      ", t, n);
			print("Pattern   ", p, m);
			print("Naive     ", m1, n);
			print("NsqrtK    ", m9, n);
			cout << "First difference at " << i << "naive said" << m1[i]
					<< " other " << m9[i] << endl;
			print("Text diff ", t + i, n - i);
			print("Pattern   ", p, m);
			exit(1);
		}
		delete[] m9;
	}

	if (false)
	{ // Randomized O(nklogn)
		if (m * k >= 500000) {
			cout << " \t -";
		} else {
			resetTime();
			int *m3 = randomizedKMismatches(t, n, p, m, k, 3);
			cout << " \t " << getTime();
			//	printTime("Randomized");
			convertTo01(m1, n, k);
			if (memcmp(m3, m1, n * sizeof(int)) != 0) {
				cout << "Difference between naive and randomized version"
						<< endl;
				exit(0);

			}
			delete[] m3;
		}
	}
	delete[] m1;
}

// randomly pick a pattern from the text
int *pick_patt(int *t, int n, int m) {
	int *p = new int[m];
	int pos = n == m ? 0 : rand() % (n - m);
	memcpy(p, t + pos, m * sizeof(int));
//	cout << "Picked pattern starting at pos " << pos << endl;
	return p;
}

void plant_and_test(int sigma, int *t, int n, int * p, int m, int *K, int nk,
		int plantingType) {
	if (plantingType == 2) {
		for (int i = 0; i + m <= n; i += m) {
			// plant pattern
			memcpy(t + i, p, m * sizeof(int));
		}
	}

	int *tcopy = new int[n];
	for (int i = 0; i < n; ++i)
		tcopy[i] = t[i];
	for (int i = 0; i < nk; ++i) {
		int k = K[i] * m / 100;
		// mutate k positions
		if (plantingType == 2) {
			for (int r = 0; r + m <= n; r += m) {
				for (int j = 0; j < k; ++j) {
					int pos = rand() % m;
					t[r + pos] = (t[r + pos] + 1) % sigma;
				}
			}
		}
		cout << sigma << " \t " << plantingType << " \t ";
		//		fflush(stdout);
		test(t, n, p, m, k);
		cout << " \t OK" << endl; //>printf("OK\n");
		memcpy(t, tcopy, n * sizeof(int));
	}

	delete[] tcopy;
}

void generateAndTestData(int sigma, int n, int m, int *K, int nk,
		int plantingType) {
	if (m >= n)
		return;
	int *t = new int[n];
	for (int i = 0; i < n; ++i)
		t[i] = rand() % sigma;
	int *p = pick_patt(t, n, m);
	plant_and_test(sigma, t, n, p, m, K, nk, plantingType);
	delete[] t;
	delete[] p;
}

void generate_tests(int sigma, int n, int m, int plant) {
	srand(123 + sigma + n + m + plant);
	generateAndTestData(sigma, n, m, K, sizeof(K) / sizeof(int), plant);
}

void generate_tests() {
	int sig[] = { 4, 20 };
	int n[] = { 1000, 10000, 100000, 1000000 };
	int m[] = { 100, 200, 1000, 2000 };
	int plant[] = { 1, 2 };

	cout
			<< "Sigma \t Plant \t n \t m \t k \t Naive \t Abrah \t NK-m2 \t NK-sa \t NK-st \t Nsqm2 \t Nsqsa \t Nsqst \t Randomized"
			<< endl;
	for (int a = 0, b = sizeof(sig) / sizeof(int); a < b; ++a)
		for (int i = 0, l = sizeof(n) / sizeof(int); i < l; ++i)
			for (int j = 0, o = sizeof(m) / sizeof(int); j < o; ++j)
				for (int ptype = 0, c = sizeof(plant) / sizeof(int); ptype < c; ++ptype)
					generate_tests(sig[a], n[i], m[j], plant[ptype]);
}

void test_conv() {
	int a[] = { 1, 2, 3 };
	int b[] = { 1, 2 };
	int n = sizeof(a) / sizeof(int);
	int m = sizeof(b) / sizeof(int);
	int k = n + m - 1;
	int c[k];
	print("A ", a, n);
	print("B ", b, m);
	Convolution conv(k);
	conv.convolve(a, n, b, m, c);
	print("C ", c, n);
}

void test_preproc() {
	int p[] = { 1, 2, 3, 3, 2, 3, 1 };
	int m = sizeof(p) / sizeof(int);
	int t[] = { 1, 2, 1, 3, 2, 3, 2 };
	int n = sizeof(t) / sizeof(int);

	print("pattern ", p, m);
	print("text    ", t, n);

	stpreproc pp(p, m);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < m; ++j)
			pp.lcp(i, j);

	for (int i = 0; i < n; ++i) {
		int suffix;
		int l = pp.longestMatch(t + i, n - i, suffix);
		cout << "Longest match for start pos " << i << " is " << l
				<< " - pattern suffix " << suffix << endl;
	}
}

unsigned char* readData(const char * const filename, int& n) {
	struct stat fileInfo;
	FILE *file;
	if (stat(filename, &fileInfo)) {
		cout << "Unable to get stat of file " << filename << endl;
		return NULL;
	}
	n = fileInfo.st_size;
	unsigned char *result = new unsigned char[n];
	if (!(file = fopen(filename, "r"))) {
		cout << "Unable to open file " << filename << endl;
		return NULL;
	}
	rewind(file);
	if (n > fread(result, sizeof(uchar), n, file)) {
		cout << "Error reading file " << filename << endl;
		fclose(file);
		delete[] result;
		return NULL;
	}
	fclose(file);
	return result;
}

void test_from_file(char *file_name) {
	int n = 0;
	cout << "reading file " << file_name << endl;
	unsigned char *c = readData(file_name, n);
	if (c == NULL)
		return;
	int *t = new int[n];
	int min = t[0];
	for (int i = 0; i < n; ++i) {
		t[i] = c[i];
		if (t[i] < min)
			min = t[i];
	}
	for (int i = 0; i < n; ++i)
		t[i] -= min;

	delete[] c;
	cout << "read " << n << " chars" << endl;

	int M[] = { 100, 200, 1000, 2000 };
	for (unsigned j = 0; j < sizeof(M) / sizeof(int); ++j) {
		int m = M[j];
		for (unsigned i = 0; i < sizeof(K) / sizeof(int); ++i)
			for (int r = 0; r < 1; ++r) {
				cout << file_name << " \t ";
				int *p = pick_patt(t, n, m);
				test(t, n, p, m, (int) (K[i] * m / 100));
				cout << " OK" << endl;
				delete[] p;
			}
	}
	delete[] t;

}

int main(int argc, char **argv) {
	//	test_preproc();
	//	return 0;

	//	test_conv();
	//	return 0;

	//	generate_tests();
	//	return 0;

	if (argc < 2) {
		cout << "Arguments: fileNames" << endl;
	}
	srand(123);
	for (int i = 1; i < argc; ++i)
		test_from_file(argv[i]);
	return 0;

	if (argc < 5) {
		cout << "Arguments: sigma n m plantType" << endl;
	} else {
		int sigma = atoi(argv[1]);
		int n = atoi(argv[2]);
		int m = atoi(argv[3]);
		int plant = atoi(argv[4]);
		generate_tests(sigma, n, m, plant);
	}
	return 0;

	if (argc < 3) {
		cout << "Arguments: input_file output_file" << endl;
		return 0;
	}
	int *t, *p;
	int m, n, k;
	readInput(argv[1], t, n, p, m, k);
	test(t, n, p, m, k);
	return 0;

	print("Text    ", t, n);
	print("Pattern ", p, m);
	cout << "k = " << k << endl;

	int *m1 = naiveCountMismatches(t, n, p, m);
	int *m2 = abrahamsonCountMismatches(t, n, p, m);
	int *m3 = randomizedKMismatches(t, n, p, m, k);
	int flag[n];
	for (int i = 0; i < n; ++i)
		flag[i] = 1;

	//		m2preproc pp(p, m);
	//		sapreproc pp(p, m);
	stpreproc pp(p, m);
	int *m4 = nkCountMismatches(t, n, p, m, k, flag, pp);

	int *m5 = nsqrtKMismatches(t, n, p, m, k, 0);
	print("Naive  ", m1, n);
	print("Abrah  ", m2, n);
	print("Rando  ", m3, n);
	print("NK     ", m4, n);
	print("NsqrtK ", m5, n);
	trim(m1, n, k + 1);
	for (int i = 0; i < n; ++i)
		if (m1[i] != m4[i]) {
			cout << "Error in alignment at position " << i << " naive said "
					<< m1[i] << " other said " << m4[i] << endl;
			exit(0);
		}
	cout << "OK" << endl;

	//	test02(n);
	//	int T[] {1, 2, 0, 1, 2, 1, 0, 1, 2, 1};
	//	int P[] {1, 0, 2};
	//	int N = sizeof(T) / sizeof(int);
	//	int M = sizeof(P) / sizeof(int);
	//	int MAX = max(N, M);
	//	int* C = new int[MAX];
	//	convolution(T, N, P, M, C);
	//	print("T ", T, N);
	//	print("P ", P, M);
	//	print("Result ", C, MAX);
	//	instance = kmis(t, n, p, m);
	//	matches = instance.find(nm);
	//	FILE *f = fopen(argv[2], "w");
	//	output(argv[2], matches, nm);
	//	fclose(f);
	delete[] t;
	delete[] p;
	delete[] m1;
	delete[] m2;
	delete[] m3;
	delete[] m4;
	delete[] m5;
	return 0;
}
