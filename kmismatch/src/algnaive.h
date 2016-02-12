#ifndef __algnaive_h__
#define __algnaive_h__ 1
#include "kutil.h"

/**
 * Counts for every position i in the text the number of mismatches between the
 * text and the pattern aligned at position i of the text. This is the naive O(n m)
 * method.
 *
 * Returns an array of the same size as the text.
 */
int *naiveCountMismatches(int *t, int n, int *p, int m) {
	int *mism = new int[n];
	memset(mism, 0, n * sizeof(int));
	for (int i = 0; i < n; ++i)
		for (int j = 0, l = min(n - i, m); j < l; ++j)
			mism[i] += t[i + j] == p[j] ? 0 : 1;
	return mism;
}

#endif
