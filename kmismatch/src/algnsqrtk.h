#ifndef __algnsqrtk_h__
#define __algnsqrtk_h__ 1
#include "kutil.h"
#include "algnk.h"

void charFreq(int *t, int n, int *alpha, int r, int *F) {
	for (int i = 0; i < n; ++i) {
		int c = binarySearch(t[i], alpha, r);
		if (c == -1)
			c = r;
		F[c]++;
	}
}

/**
 * Counts for every position i in the text, the number of mismatches between the
 * text and the pattern aligned at position i, if the number is less than k, otherwise,
 * a number greater than k is reported. This is a simplified version of Amir's algorithm.
 *
 * Returns an array of the same size as the text.
 */
int *nsqrtKMismatches(int *t, int n, int *p, int m, int k, int preprocessorType) {
	// alphabet and frequency of each character in the pattern
	int alpha[m];
	int f[m];
	int r = alphaFreq(p, m, alpha, f);

	// compute frequency of characters in the text
	int F[r + 1];
	memset(F, 0, sizeof(int) * (r + 1));
	charFreq(t, n, alpha, r, F);

	// sort characters by increasing frequency in the text!
	mypair a[r];
	for (int i = 0; i < r; ++i) {
		a[i].a = i;
		a[i].b = F[i];
	}
	qsort(a, r, sizeof(mypair), pair_b_cmp);

	// try to fill a knapsack of size 2k
	int knapsack = 2 * k;
	int relevantInst[r];
	memset(relevantInst, 0, sizeof(int) * r);
	long maxCost = round(5 * (long) n * (double) sqrt(k * intLog(k + 1)));
	int d = 0;
	long cost = 0;
	int items = 0;
	for (; cost < maxCost && items < knapsack && d < r; ++d) {
		int index = a[d].a;
		int toUse = min(f[index], knapsack - items);
		toUse = min(toUse, 1 + (maxCost - cost) / F[index]);
		cost += (long)F[index] * toUse;
		items += toUse;
		relevantInst[index] = toUse;
	}

	bool knapsackFilled = items == knapsack;
	if (!knapsackFilled) {
		int last = a[d - 1].a;
		if (relevantInst[last] < f[last]) {
			// if last item is only partially included in knapsack,
			// we leave it to be counted by convolution
			relevantInst[last] = 0;
			--d;
		}
	}

	// compute offsets and use counting for items in knapsack
	int *relevantOffset[r];
	for (int i = 0; i < r; ++i)
		relevantOffset[i] = new int[knapsack];
	int found[r];
	memset(found, 0, sizeof(int) * r);
	for (int i = 0; i < m; ++i) {
		int index = binarySearch(p[i], alpha, r);
		if (found[index] < relevantInst[index])
			relevantOffset[index][found[index]++] = i;
	}
	int *matches = new int[n];
	memset(matches, 0, n * sizeof(int));
	count(alpha, r, t, n, relevantInst, relevantOffset, matches);
	for (int i = 0; i < r; ++i)
		delete[] relevantOffset[i];

	if (knapsackFilled) {
		// for each position with more than k counts, check alignment using a O(k)/position algorithm
		int *flag = matches;
		for (int i = 0; i < n; ++i)
			if (matches[i] < k)
				flag[i] = 0;
			else
				flag[i] = 1;
		for (int i = n - m + 1; i < n; ++i)
			flag[i] = 1;
		preprocessor pp;
		if (preprocessorType == 0) {
			pp = m2preproc(p, m);
		} else if (preprocessorType == 1) {
			pp = sapreproc(p, m);
		} else {
			pp = stpreproc(p, m);
		}
		int *mismatches = nkCountMismatches(t, n, p, m, k, flag, pp);
		for (int i = 0; i < n; ++i)
			if (!flag[i])
				mismatches[i] = k + 1;
		delete[] matches;
		return mismatches;
	} else {
		// count the remaining characters using convolutions
		int remaining = r - d;
		int rest[remaining];
		for (int i = d; i < r; ++i) {
			rest[i - d] = alpha[a[i].a];
		}
		convolutions_fftw(rest, remaining, t, n, p, m, matches);
		//		convolutions(rest, remaining, t, n, p, m, matches);
		int *mismatches = matches;
		for (int i = 0; i < n; ++i) {
			int mis = min(m, n - i) - matches[i];
			mismatches[i] = mis;
		}
		return mismatches;
	}
}
#endif
