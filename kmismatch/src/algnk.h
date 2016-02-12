#ifndef __algnk_h__
#define __algnk_h__ 1
#include "kutil.h"
#include "SSTree.h"
#include "radix.h"

struct entry {
	unsigned hCode;
	int suffix, length;
	entry(unsigned hCode, int suffix, int length) :
		hCode(hCode), suffix(suffix), length(length) {

	}
};

class preprocessor {
public:
	virtual int longestMatch(int *t, int n, int &suffix) {
		return 0;
	}
	virtual int lcp(int s1, int s2) {
		return 0;
	}
};

class nopreproc: public preprocessor {
private:
	int *p;
	int m;

public:
	nopreproc(int *p, int m) :
		p(p), m(m) {

	}
	int longestMatch(int *t, int n, int &suffix) {

		// naive
		int max = 0;
		for (int i = 0; i < m; ++i) {
			int *c = t;
			int j = 0;
			while (i + j < m && j < n && *c == p[i + j]) {
				++j;
				++c;
			}
			if (j > max) {
				max = j;
				suffix = i;
			}
		}
		return max;
	}

	int lcp(int s1, int s2) {
		//naive
		int l = 0;
		while (s1 + l < m && s2 + l < m && p[s1 + l] == p[s2 + l])
			++l;
		return l;
	}

};

class m2preproc: public preprocessor {
private:
	int *p;
	int m;
	int **a;
	vector<list<entry> > hStart;
	int hSize;

	int **buildLCPMatrix(int *p, int m) {
		int **a = new int*[m + 1];
		for (int i = 0; i <= m; ++i)
			a[i] = new int[m + 1];
		for (int i = 0; i <= m; ++i)
			a[i][m] = a[m][i] = 0;
		for (int i = m - 1; i >= 0; --i)
			for (int j = m - 1; j >= 0; --j)
				if (p[i] == p[j])
					a[i][j] = 1 + a[i + 1][j + 1];
				else
					a[i][j] = 0;
		return a;
	}

	unsigned updateHash(unsigned hc, int ch) {
		return hc * 31 + ch;
	}

	bool isPrime(int n) {
		for (int i = 2; i * i <= n; ++i)
			if (n % i == 0)
				return false;
		return true;
	}

	void buildHashMap(int *p, int m) {
		hSize = 2 * m * m + 1;
		while (!isPrime(hSize))
			++hSize;
		hStart = vector<list<entry> > (hSize);
		for (int i = 0; i < m; ++i) {
			unsigned int hc = hSize;
			for (int j = 0, e = m - i; j < e; ++j) {
				hc = updateHash(hc, p[i + j]);
				int e = hc % hSize;
				hStart[e].push_back(entry(hc, i, j + 1));
			}
		}
	}

	bool match(int *p, int m, int start, int *t, int n, int len) {
		if (start + len > m || n < len)
			return false;
		p += start;
		while (len && *p++ == *t++)
			--len;
		return len == 0;
	}

public:
	m2preproc(int *p, int m) :
		p(p), m(m) {
		a = buildLCPMatrix(p, m);
		buildHashMap(p, m);
	}

	~m2preproc() {
		for (int i = 0; i <= m; ++i)
			delete[] a[i];
		delete[] a;
	}

	int longestMatch(int *t, int n, int &suffix) {
		/*/////////////////////////
		 // naive
		 int max = 0;
		 for (int i = 0; i < m; ++i) {
		 int *c = t;
		 int j = 0;
		 while (i + j < m && j < n && *c == p[i + j]) {
		 ++j;
		 ++c;
		 }
		 if (j > max) {
		 max = j;
		 suffix = i;
		 }
		 }
		 *//////////////////////////

		//		return max;
		unsigned hc = hSize;
		int i = 0;
		while (i < n && i < m) {
			unsigned nextHc = updateHash(hc, t[i]);
			bool stop = true;
			list<entry>& l = hStart[nextHc % hSize];
			for (list<entry>::iterator it = l.begin(); it != l.end(); ++it) {
				if (it->hCode == nextHc && it->length == i + 1) {
					stop = false;
					break;
				}
			}
			if (stop)
				break;
			hc = nextHc;
			++i;
		}
		if (i == 0)
			return 0;
		do {
			list<entry>& a = hStart[hc % hSize];
			for (list<entry>::iterator it = a.begin(); it != a.end(); ++it) {
				int start = it->suffix;
				if (it->hCode == hc && match(p, m, start, t, n, i)) {
					suffix = start;
					return i;
				}
			}
			if (i == 1) {
				return 0;
			} else {
				--i;
				hc = hSize;
				for (int j = 0; j < i; ++j)
					hc = updateHash(hc, t[i]);
			}
		} while (1);

	}

	/**
	 * Returns the length of the longest common prefix of two suffixes
	 */
	int lcp(int s1, int s2) {
		//		naive
		//		int l = 0;
		//		while (s1 + l < m && s2 + l < m && p[s1 + l] == p[s2 + l])
		//			++l;
		if (s1 < 0 || s2 < 0 || s1 > m || s2 > m) {
			cout << "AUCH " << s1 << " " << s2 << endl;
			exit(0);
		}
		int ans = a[s1][s2];
		//		if (ans != l)
		//			cout << endl << "BUBA - LCP should be " << l << " not " << ans
		//					<< endl;
		return ans;
	}

};
/**
 * Counts for every position i in the text for which flag[i]=1, the number of mismatches between the
 * text and the pattern aligned at position i. If the number is less than k, it is reported, if it's
 * greater than k, then k+1 is reported.
 *
 * Returns an array of the same size as the text.
 */
int *nkCountMismatches(int *t, int n, int *p, int m, int k, int *flag,
		preprocessor &pp) {
	int *mism = new int[n];
	memset(mism, 0, n * sizeof(int));
	list<int> q;
	for (int i = 0; i < n; ++i) {
		int sm = 0;
		int l = pp.longestMatch(t + i, n - i, sm);
		for (int j = 0; j <= l; ++j) {
			if (i + j < n && flag[i + j])
				q.push_back(i + j);
		}
		for (list<int>::iterator it = q.begin(); it != q.end();) {
			int alignStart = i - *it;
			if (alignStart < m) {
				int j = alignStart < 0 ? -alignStart : 0;
				do {
					j += pp.lcp(alignStart + j, sm + j);
					if (j < l) {
						if (alignStart + j < m)
							++mism[*it];
					} else if (alignStart + l < m && i + l < n && t[i + l]
							!= p[alignStart + l])
						mism[*it]++;

					if (mism[*it] > k)
						break;
					++j;
				} while (j <= l && alignStart + j < m);
			}
			if (mism[*it] > k || alignStart >= m)
				it = q.erase(it);
			else
				it++;
		}
		i += l;
	}
	return mism;
}

class sapreproc: public preprocessor {
private:
	int *p;
	int m;
	unsigned *sa;
	int *pos;
	int *lcpTable;
	int **minRange;

	int charAt(int pos) {
		return pos < m ? p[pos] : -1;
	}

	void narrowInterval(int lo, int hi, int key, int depth, int &nlo, int &nhi) {
		while (lo < hi) {
			int k = (lo + hi) >> 1;
			int c = charAt(sa[k] + depth);
			if (c < key) {
				lo = k + 1;
			} else if (c > key) {
				hi = k;
			} else {
				if (charAt(sa[lo] + depth) == charAt(sa[hi - 1] + depth))
					break;
				int tmp;
				narrowInterval(lo, k, key, depth, nlo, tmp);
				narrowInterval(k + 1, hi, key, depth, tmp, nhi);
				return;
			}
		}
		nlo = lo;
		nhi = hi;
	}

	int rmq(int l, int r) {
		if (r < l) {
			int tmp = l;
			l = r;
			r = tmp;
		}
		int bit = intLog(r - l);
		return min(minRange[l][bit], minRange[r - (1 << bit)][bit]);
	}

	void buildSA() {
		uchar pc[m];
		for (int i = 0; i < m; ++i) {
			pc[i] = p[i];
		}
		Radix r(pc, m);
		sa = r.build();
	}

	void buildLcpTable() {
		lcpTable = new int[m];
		int h = 0;
		for (int i = 0; i < m; ++i) {
			int k = pos[i] + 1;
			if (k < m) {
				int j = sa[k];
				int gap = min(m - i, m - j);
				while (h < gap && p[i + h] == p[j + h])
					++h;
				lcpTable[pos[i]] = h;
			} else {
				lcpTable[pos[i]] = 0;
			}
			if (h > 0)
				--h;
		}
	}

	void preprocessMinRange() {
		minRange = new int*[m];
		int size = 2 + intLog(m);
		for (int i = m - 1; i >= 0; --i) {
			minRange[i] = new int[size];
			minRange[i][0] = lcpTable[i];
			for (int b = 0, eb = 1; i + eb < m; ++b, eb <<= 1)
				minRange[i][b + 1] = min(minRange[i][b], minRange[i + eb][b]);
		}
	}

	void init() {
		buildSA();
		pos = new int[m];
		for (int i = 0; i < m; ++i)
			pos[sa[i]] = i;
		buildLcpTable();
		preprocessMinRange();
	}

	void clear() {
		delete[] sa;
		delete[] pos;
		delete[] lcpTable;
		for (int i = 0; i < m; ++i)
			delete[] minRange[i];
		delete[] minRange;
	}

public:
	sapreproc(int *p, int m) :
		p(p), m(m) {
		init();
	}

	~sapreproc() {
		clear();
	}
	int longestMatch(int *t, int n, int &suffix) {
/*
		 // naive
		 int max = 0;
		 for (int i = 0; i < m; ++i) {
		 int *c = t;
		 int j = 0;
		 while (i + j < m && j < n && *c == p[i + j]) {
		 ++j;
		 ++c;
		 }
		 if (j > max) {
		 max = j;
		 suffix = i;
		 }
		 }
*/
		int i = 0;
		int lo = 0;
		int hi = m;
		do {
			int nlo, nhi;
			narrowInterval(lo, hi, t[i], i, nlo, nhi);
			if (nlo < nhi) {
				++i;
				lo = nlo;
				hi = nhi;
			} else
				break;
		} while (i < n);
		suffix = sa[lo];
//				if (i != max)
//					printf("Problem - longest match should be %d but is %d\n", max, i);
		return i;
		//		return max;
	}

	int lcp(int s1, int s2) {
		//naive
//		 int l = 0;
//		 while (s1 + l < m && s2 + l < m && p[s1 + l] == p[s2 + l])
//		 ++l;

		int r = m - s1;
		if (s1 == m || s2 == m)
			r = 0;
		else if (s1 != s2)
			r = rmq(pos[s1], pos[s2]);
//				if (r != l)
//					printf("Problem - lcp should be %d but is %d\n", l, r);
		return r;
	}

};

class stpreproc: public preprocessor {
private:
	int *p;
	unsigned char *pc;
	int m;
	SSTree *sst;

public:
	stpreproc(int *p, int m) :
		p(p), m(m) {
		pc = new unsigned char[m + 1];
		for (int i = 0; i < m; ++i)
			pc[i] = (unsigned char) ('a' + p[i]);
		pc[m] = '\0';
		sst = new SSTree(pc, (unsigned long) (m + 1));
		//		sst->PrintTree(sst->root(), 0);
	}

	~stpreproc() {
		delete sst;
		delete[] pc;
	}

	int longestMatch(int *t, int n, int &suffix) {
		/*
		 // naive
		 int max = 0;
		 int suffix2 = 0;
		 for (int i = 0; i < m; ++i) {
		 int *c = t;
		 int j = 0;
		 while (i + j < m && j < n && *c == p[i + j]) {
		 ++j;
		 ++c;
		 }
		 if (j > max) {
		 max = j;
		 suffix2 = i;
		 }
		 }
		 */
		ulong node = sst->root();
		ulong child = 0;
		int i = 0;
		while (i < n && 0 != (child = sst->child(node, 'a' + t[i]))) {
			unsigned long ps, len;
			sst->edge(child, ps, len);
			int *label = p + ps;
			unsigned int j = 0;
			unsigned int lim = m - ps;
			for (; j < len && j < lim && i < n && (t[i] == label[j]); ++j, ++i)
				;
			node = child;
			if (j < len)
				break;
		}
		//		cout << "Path label " << sst->pathlabel(node) << endl;
		suffix = sst->textpos(node);
		/*
		 if (i != max) {
		 printf("\nProblem - longest match should be %d but is %d\n", max, i);
		 printf("Suffix should be %d and is %d\n", suffix2, suffix);
		 print("text: ", t, n);
		 print("patt: ", p + suffix, m - suffix);
		 }*/
		return i;
		//		return max;
	}

	int lcp(int s1, int s2) {
		//naive
		//		int l = 0;
		//		while (s1 + l < m && s2 + l < m && p[s1 + l] == p[s2 + l])
		//			++l;

		int r = s1 == s2 ? m - s1 : sst->lce(s1, s2);
		/*if (r != l) {
		 printf("\nProblem - lcp of %d and %d should be %d but is %d\n", s1,
		 s2, l, r);
		 print("first :", p + s1, m - s1);
		 print("second:", p + s2, m - s2);
		 exit(0);
		 }*/
		return r;
	}

};

#endif
