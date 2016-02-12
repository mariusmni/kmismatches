#ifndef __kutil_h__
#define __kutil_h__ 1

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <vector>
#include <string>
#include <cmath>
//#include "convolution.h"
#define KDEBUG 0
#include <fftw3.h>

using namespace std;
//using namespace fftwpp;

int int_cmp(const void *a, const void *b) {
	return *((int *) a) - *((int *) b);
}

/**
 * Computes the alphabet of a text and the frequency of each character in the text.
 * Returns the size of the alphabet. Returns a sorted list of distinct characters in alpha
 * and their frequencies in freq.
 */
int alphaFreq(int *p, int m, int *alpha, int *freq) {
	memcpy(alpha, p, m * sizeof(int));
	qsort(alpha, m, sizeof(int), int_cmp);

	memset(freq, 0, m * sizeof(int));
	int j = -1;
	int r = 0;
	for (int i = 0; i < m; ++i)
		if (alpha[i] == j)
			freq[r - 1]++;
		else {
			freq[r] = 1;
			alpha[r] = alpha[i];
			j = alpha[i];
			r++;
		}
	return r;
}

void print(string msg, int *a, int n) {
	cout << msg;
	for (int i = 0; i < n; ++i)
		cout << a[i] << " ";
	cout << endl;
}

/*
void print(string msg, Complex *a, int n) {
	cout << msg;
	for (int i = 0; i < n; ++i)
		cout << a[i] << " ";
	cout << endl;
}
*/

/*
void convolution(Complex *f, Complex *g, ImplicitConvolution &C, int M, int *c,
		int n) {
	//	if (DEBUG) {
	//		print("F ", f, M);
	//		print("G ", g, M);
	//	}
	C.convolve(f, g);
	//	if (DEBUG) {
	//		print("F*G ", f, M);
	//		print("G after ", g, M);
	//	}
	for (int i = 0; i < n; ++i)
		c[i] = round(f[M - i - 1].re);
	//	if (DEBUG)
	//		print("ci = ", c, n);
}

void convolution(int *a, int n, int *b, int m, int *c, Complex *f, Complex *g,
		ImplicitConvolution &C, int M) {
	do {
		for (int i = 0; i < M; ++i)
			f[i] = g[i] = Complex(0, 0);
		for (int i = 0; i < m; ++i)
			g[i] = Complex(b[i], 0);
		int mm = min(M, n);
		for (int i = max(M - n, 0), j = mm - 1; i < M; ++i, --j)
			f[i] = Complex(a[j], 0);
		convolution(f, g, C, M, c, mm);

		n -= m;
		a += m;
		c += m;
	} while (n > m);
}*/

/**
 * Computes the convolution of vectors a of size n
 * and b of size m and returns the result in vector c.
 * c should have a size at least equal to max(m,n)
 *
 * The convolution is defined as:
 * c[i] = sum_{j=0,m-1}(a[i+j] * b[j])
 */
/*
void convolution(int *a, int n, int *b, int m, int *c) {
	if (KDEBUG) {
		print("Convolution of ", a, n);
		print("and ", b, m);
	}
	int M = (n > 2 * m) ? 2 * m : max(m, n);
	Complex *f = ComplexAlign(M);
	Complex *g = ComplexAlign(M);
	ImplicitConvolution C(M);
	//	ExplicitConvolution EC()

	//	int *originalC = c;
	//	int originalN = n;
	convolution(a, n, b, m, c, f, g, C, M);

	//	if (KDEBUG)
	//		print("c = ", originalC, originalN);
	deleteAlign(g);
	deleteAlign(f);
}*/

int binarySearch(int key, int *a, int n) {
	int i = 0, j = n - 1;
	while (i <= j) {
		int k = (i + j) >> 1;
		if (a[k] == key)
			return k;
		if (a[k] < key)
			i = k + 1;
		else
			j = k - 1;
	}
	return -1;
}

struct mypair {
	int a, b;
};

int pair_b_rcmp(const void *a, const void *b) {
	return ((mypair *) b)->b - ((mypair *) a)->b;
}

int pair_b_cmp(const void *a, const void *b) {
	return ((mypair *) a)->b - ((mypair *) b)->b;
}

//int min(int a, int b) {
//	return a < b ? a : b;
//}

void printComplex(char *msg, fftw_complex *a, int n) {
	cout << msg;
	for (int i = 0; i < n; ++i)
		cout << "(" << a[i][0] << ", " << a[i][1] << ") ";
	cout << endl;
}

class Convolution {
private:
	int N;
	fftw_complex *fw;
	fftw_complex *bw;
	fftw_complex *t;
	fftw_plan forward;
	fftw_plan backward;
public:
	Convolution(int N) :
		N(N) {
		fw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		bw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
		forward = fftw_plan_dft_1d(N, fw, fw, FFTW_FORWARD, FFTW_ESTIMATE);
		backward = fftw_plan_dft_1d(N, bw, bw, FFTW_BACKWARD, FFTW_ESTIMATE);
	}
	~Convolution() {
		fftw_destroy_plan(forward);
		fftw_destroy_plan(backward);
		fftw_free(fw);
		fftw_free(bw);
		fftw_free(t);
	}

	void mult(fftw_complex &x, fftw_complex &y, fftw_complex &z) {
		double re = x[0] * y[0] - x[1] * y[1];
		double im = x[0] * y[1] + x[1] * y[0];
		z[0] = re;
		z[1] = im;
	}

	void transform(int *b, int m) {
		// transform b
		for (int i = 0, j = m - 1; i < m; ++i, --j) {
			fw[i][0] = b[j];
			fw[i][1] = 0;
		}
		for (int i = m; i < N; ++i)
			fw[i][0] = fw[i][1] = 0;
		fftw_execute(forward);
		for (int i = 0; i < N; ++i) {
			t[i][0] = fw[i][0];
			t[i][1] = fw[i][1];
		}
	}

	void convolveWithTransformed(int *a, int n, int m, int *c) {
		// transform a
		for (int i = 0; i < n; ++i) {
			fw[i][0] = a[i];
			fw[i][1] = 0;
		}
		for (int i = n; i < N; ++i)
			fw[i][0] = fw[i][1] = 0;
		fftw_execute(forward);

		// multiply
		for (int i = 0; i < N; ++i)
			mult(t[i], fw[i], bw[i]);

		// reverse transform
		fftw_execute(backward);

		for (int i = 0; i < n; ++i)
			c[i] = round(bw[i + m - 1][0] / N);
	}

	void convolve(int *a, int n, int *b, int m, int *c) {
		// transform a
		for (int i = 0; i < n; ++i) {
			fw[i][0] = a[i];
			fw[i][1] = 0;
		}
		for (int i = n; i < N; ++i)
			fw[i][0] = fw[i][1] = 0;
		//		printComplex("a            ", fw, N);

		fftw_execute(forward);
		for (int i = 0; i < N; ++i) {
			bw[i][0] = fw[i][0];
			bw[i][1] = fw[i][1];
		}

		// transform b
		for (int i = 0, j = m - 1; i < m; ++i, --j) {
			fw[i][0] = b[j];
			fw[i][1] = 0;
		}
		for (int i = m; i < N; ++i)
			fw[i][0] = fw[i][1] = 0;
		//		printComplex("b            ", fw, N);
		fftw_execute(forward);
		//		printComplex("a'           ", bw, N);
		//		printComplex("b'           ", fw, N);

		// multiply
		for (int i = 0; i < N; ++i)
			mult(fw[i], bw[i], bw[i]);
		//		printComplex("a' * b'      ", bw, N);

		// reverse transform
		fftw_execute(backward);
		//		printComplex("(a' * b')^-1 ", bw, N);

		for (int i = 0; i < n; ++i)
			c[i] = round(bw[i + m - 1][0] / N);
	}

};

void convolution_fftw(int *a, int n, int *b, int m, int *c, int M,
		Convolution &conv) {
	conv.transform(b, m);
	for (int i = 0; i < n; i += m) {
		int size = min(n - i, M);
		conv.convolveWithTransformed(a + i, size, m, c + i);
	}
}

/**
 * For every character in chr computes the number of matches between t and p
 * using convolution, and adds that number to the relevant position in matches
 */
void convolutions_fftw(int *chr, int size, int *t, int n, int *p, int m,
		int *matches) {
	int P[m];
	int M = (n > 2 * m) ? 2 * m : max(m, n);
	int T[M];
	int C[M];

	Convolution conv(M + m - 1);

	for (int a = 0; a < size; ++a) {
		int c = chr[a];
		for (int j = 0; j < m; ++j)
			P[j] = (p[j] == c ? 1 : 0);
		conv.transform(P, m);
		for (int i = 0; i < n; i += m) {
			int size = min(n - i, M);
			for (int j = 0; j < size; ++j)
				T[j] = (t[i + j] == c ? 1 : 0);
			conv.convolveWithTransformed(T, size, m, C);
			for (int j = 0, l = min(size, m); j < l; ++j)
				matches[i + j] += C[j];
		}
	}
}

/**
 * For every character in chr computes the number of matches between t and p
 * using convolution, and adds that number to the relevant position in matches
 */
/*
void convolutions(int *chr, int size, int *t, int n, int *p, int m,
		int *matches) {
	int T[n];
	int C[n];
	int P[m];
	int M = (n > 2 * m) ? 2 * m : max(m, n);

	Complex *f = ComplexAlign(M);
	Complex *g = ComplexAlign(M);
	ImplicitConvolution iC(M);

	for (int i = 0; i < size; ++i) {
		int c = chr[i];
		for (int j = 0; j < n; ++j)
			T[j] = (t[j] == c ? 1 : 0);
		for (int j = 0; j < m; ++j)
			P[j] = (p[j] == c ? 1 : 0);
		convolution(T, n, P, m, C, f, g, iC, M);
		//		convolution(T, n, P, m, C);
		for (int j = 0; j < n; ++j)
			matches[j] += C[j];
	}
	deleteAlign(g);
	deleteAlign(f);
}
*/

/**
 * For every character in alphabet counts the number of matches at every alignment
 * in t, using the counting method. The relevant number of offsets and the offsets
 * themselves, for characters in alphabet are given in offsets and offset respectively.
 * The counts are added to the corresponding positions in matches.
 */
void count(int *alphabet, int alphabetSize, int *t, int n, int *offsets,
		int **offset, int *matches) {
	for (int i = 0; i < n; ++i) {
		int c = binarySearch(t[i], alphabet, alphabetSize);
		if (c >= 0)
			for (int j = offsets[c] - 1; j >= 0; --j) {
				int align = i - offset[c][j];
				if (align >= 0)
					matches[align]++;
			}
	}
}

int intLog(int n) {
	int b = 0;
	if (n >= 0x10000) {
		b += 16;
		n >>= 16;
	}
	if (n >= 0x100) {
		b += 8;
		n >>= 8;
	}
	if (n >= 0x10) {
		b += 4;
		n >>= 4;
	}
	if (n >= 0x4) {
		b += 2;
		n >>= 2;
	}
	if (n >= 0x2) {
		b += 1;
		n >>= 1;
	}
	return b;
}
#endif
