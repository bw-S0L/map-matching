
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<algorithm>
#include<queue>
#include<cstring>
#define N 205
#define ROADSIZE 125000
#define pi 3.1415926
#define SIZE 20
using namespace std;


const double  minw = 30.640555, minj = 121.003750, wlen = 0.837223, jlen = 1.089605;
const double  xstep = wlen / (N - 5), ystep = jlen / (N - 5), fiftymeter = 0.00054;
const double sigma = 6.57, ONEW = 111000 / 2 * 1.732, lambda = 1.5, Beita = 0.01;
const double Bxs = 1.0 / (sqrt(2 * pi) * sigma);
const double ln2 = 0.69314718055995;
/*
const double Beita[31] = { 0,0.49037673,0.82918373,1.24364564,1.67079581,2.00719298,2.42513007,2.81248831,
3.15745473,3.52645392,4.09511775,4.67319795,5.41088180,6.47666590,6.29010734,7.80752112,8.09074504 ,8.08550528,
9.09405065 ,11.09090603 ,11.87752824 ,12.55107715 ,15.82820829,17.69496773,18.07655652,19.63438911,
25.40832185,23.76001877,28.43289797,32.21683062,34.56991141};*/
double Fxs = 1.0 / Beita;

struct point {
	double x = 0, y = 0;
	point() {}
	point(double w, double j) :x(w), y(j) {}
};
struct pointT {
	int t;
	double x = 0, y = 0;
	pointT() {}
	pointT(int t, double w, double j) :t(t), x(w), y(j) {}
};
struct road {
	int id, idbegin, idend, c, wtype;
	double Dx, Dy, len;
	string way;
	vector<point> p;
};
struct condidate {
	int idfi, idse;
	double Dx, Dy, B;
	condidate() {}
	condidate(pair<int, int> t, double x, double y, double b) :idfi(t.first), idse(t.second), Dx(x), Dy(y), B(b) {}
};
struct ROADTYPE {
	struct p {
		int fi = 0, se = 0, te = 0;
	};
	p r[8];
	ROADTYPE() {
		for (int i = 1; i <= 7; i++)
			r[i].fi = i;
		r[1].se = r[1].te = 5;
		r[2].se = 30;
		r[2].te = 50;
		r[3].se = 20;
		r[3].te = 50;
		r[4].se = 20;
		r[4].te = 50;
		r[5].se = 40;
		r[5].te = 80;
		r[6].se = 60;
		r[6].te = 120;
	}
	int isroad(int k, double v) {
		if (v == 0)
			return 1;
		if (v <= r[k].te && v >= r[k].se)
			return 1;
		else if (v > r[k].te) {
			if (1.0 * (v - r[k].te) / r[k].te <= 0.3)
				return 1;
		}
		else if (v < r[k].se) {
			if (1.0 * (r[k].se - v) / r[k].se <= 0.3)
				return 1;
		}
		if (k == 7 && v > r[7].te)
			return 1;
		if (k == 1 || k == 2 && v < r[k].se)
			return 1;
		return 0;
	}
};
struct PRE {
	int qian[SIZE] = { 0 };
	double F[SIZE] = { 0 };
	double v[SIZE] = { 0 };
	PRE() {}
};
struct idbegend {
	int id, dz;
	idbegend() {}
	idbegend(int id, int dz) :id(id), dz(dz) {}
};


int n, m;
road a[ROADSIZE];
vector<pair<int, int> > matrix[N][N];
vector<vector<int>>ANS;
ROADTYPE roadtype;
vector<idbegend> grid[ROADSIZE];
double D[ROADSIZE];
int vis[ROADSIZE];
double maxF = 0, maxB = 0, minF = 1000, minB = 1000;


/*****************/

/*******************///测试参数区

bool cmp(condidate x, condidate y) {
	return x.B > y.B;
}
inline int getx(double x) {
	return  (N - 5) - (x - minw) / xstep;
}
inline int gety(double y) {
	return  (y - minj) / ystep;  //(0,0)=>(minw+wlen,minj)
}
inline void getleft_up(int x, int y, double& w, double& j) {
	w = -x * xstep + minw + wlen;
	j = y * ystep + minj;
}
inline void getleft_down(int x, int y, double& w, double& j) {
	w = -(x + 1) * xstep + minw + wlen;
	j = y * ystep + minj;
}
inline void getright_up(int x, int y, double& w, double& j) {
	w = -x * xstep + minw + wlen;
	j = (y + 1) * ystep + minj;
}
inline void getright_down(int x, int y, double& w, double& j) {
	w = -(x + 1) * xstep + minw + wlen;
	j = (y + 1) * ystep + minj;
}
inline void DJS(condidate u, condidate v, double& dis, double& F,int time) {
	int i = u.idfi, j = v.idfi;
	int l = u.idse, r = v.idse;
	int t = max(a[i].wtype, a[j].wtype);
	double maxdis = 1.0 * roadtype.r[t].te/18*5*time/ONEW;
	int tag = 0;
	priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>> >ttg;
	memset(vis, 0, sizeof(vis));
	for (int k = 0; k < ROADSIZE - 1; k++)
		D[k] = 10;
	ttg.push(make_pair(0, i));
	D[i] = 0;
	int sum = 0;
	dis += sqrt((u.Dx - a[i].p[l].x) * (u.Dx - a[i].p[l].x) + (u.Dy - a[i].p[l].y) * (u.Dy - a[i].p[l].y));
	dis += sqrt((v.Dx - a[j].p[r - 1].x) * (v.Dx - a[j].p[r - 1].x) + (v.Dy - a[j].p[r - 1].y) * (v.Dy - a[j].p[r - 1].y));
	for (int k = l + 1; k < a[i].p.size(); k++)
		dis += sqrt((a[i].p[k].x - a[i].p[k - 1].x) * (a[i].p[k].x - a[i].p[k - 1].x) + (a[i].p[k].y - a[i].p[k - 1].y) * (a[i].p[k].y - a[i].p[k - 1].y));
	for (int k = 1; k < r; k++)
		dis += sqrt((a[j].p[k].x - a[j].p[k - 1].x) * (a[j].p[k].x - a[j].p[k - 1].x) + (a[j].p[k].y - a[j].p[k - 1].y) * (a[j].p[k].y - a[j].p[k - 1].y));
	while (ttg.size()) {
		pair<double, int>tmp = ttg.top();
		ttg.pop();
		if (vis[tmp.second])
			continue;
		vis[tmp.second] = 1;
		sum++;
		if (D[tmp.second] > maxdis-dis)
			break;
		if (vis[j])
			break;
		for (int k = 0; k < grid[a[tmp.second].idend].size(); k++) {
			if (grid[a[tmp.second].idend][k].dz == 1 && !vis[grid[a[tmp.second].idend][k].id]) {
				double tmpdis = 0;
				int now = grid[a[tmp.second].idend][k].id;
				if (tag) {
					for (int t = 1; t < a[tmp.second].p.size(); t++)
						tmpdis += sqrt((a[tmp.second].p[t].x - a[tmp.second].p[t - 1].x) * (a[tmp.second].p[t].x - a[tmp.second].p[t - 1].x) + (a[tmp.second].p[t].y - a[tmp.second].p[t - 1].y) * (a[tmp.second].p[t].y - a[tmp.second].p[t - 1].y));
				}
				if (D[now] > D[tmp.second] + tmpdis) {
					D[now] = D[tmp.second] + tmpdis;
					ttg.push(make_pair(D[now], now));
				}
			}
		}
		tag = 1;
	}
	if (!vis[j]) {
		dis = 1;
		return;
	}
	
	dis += D[j];
}
inline void PDJC(condidate u, condidate v, double& F, double& dis,int time) {
	int i = u.idfi, j = v.idfi;
	int l = u.idse, r = v.idse;
	if (abs(a[i].wtype - a[j].wtype) <= 1)
		F *= 2;
	F = 1;
	dis = 0;
	if (i == j) {
		if (l > r) {
			dis = 1;
			return;
		}
		if (l == r) {
			dis = sqrt((u.Dx - v.Dx) * (u.Dx - v.Dx) + (u.Dy - v.Dy) * (u.Dy - v.Dy));
		}
		else if (l < r) {
			dis += sqrt((u.Dx - a[i].p[l].x) * (u.Dx - a[i].p[l].x) + (u.Dy - a[i].p[l].y) * (u.Dy - a[i].p[l].y));
			dis += sqrt((v.Dx - a[i].p[r - 1].x) * (v.Dx - a[i].p[r - 1].x) + (v.Dy - a[i].p[r - 1].y) * (v.Dy - a[i].p[r - 1].y));
			for (int k = l + 1; k < r; k++) {
				dis += sqrt((a[i].p[k].x - a[i].p[k - 1].x) * (a[i].p[k].x - a[i].p[k - 1].x) + (a[i].p[k].y - a[i].p[k - 1].y) * (a[i].p[k].y - a[i].p[k - 1].y));
			}
		}
		return;
	}
	if (a[i].idend == a[j].idbegin) {
		dis += sqrt((u.Dx - a[i].p[l].x) * (u.Dx - a[i].p[l].x) + (u.Dy - a[i].p[l].y) * (u.Dy - a[i].p[l].y));
		dis += sqrt((v.Dx - a[j].p[r - 1].x) * (v.Dx - a[j].p[r - 1].x) + (v.Dy - a[j].p[r - 1].y) * (v.Dy - a[j].p[r - 1].y));
		for (int k = l + 1; k < a[i].p.size(); k++)
			dis += sqrt((a[i].p[k].x - a[i].p[k - 1].x) * (a[i].p[k].x - a[i].p[k - 1].x) + (a[i].p[k].y - a[i].p[k - 1].y) * (a[i].p[k].y - a[i].p[k - 1].y));
		for (int k = 1; k < r; k++)
			dis += sqrt((a[j].p[k].x - a[j].p[k - 1].x) * (a[j].p[k].x - a[j].p[k - 1].x) + (a[j].p[k].y - a[j].p[k - 1].y) * (a[j].p[k].y - a[j].p[k - 1].y));
	}
	else {
		DJS(u, v, dis, F,time);
	}
}

inline void addroad(int i) {
	for (int j = 1; j < a[i].p.size(); j++) {
		double x0 = a[i].p[j - 1].x, y0 = a[i].p[j - 1].y, x1 = a[i].p[j].x, y1 = a[i].p[j].y;
		if (getx(x0) == getx(x1)) {
			int x = getx(x0);
			int y = min(gety(y0), gety(y1)), ym = max(gety(y0), gety(y1));
			for (; y <= ym; y++)
				matrix[x][y].push_back(make_pair(a[i].id, j));
			continue;
		}
		if (gety(y0) == gety(y1)) {
			int y = gety(y0);
			int x = min(getx(x0), getx(x1)), xm = max(getx(x0), getx(x1));
			for (; x <= xm; x++)
				matrix[x][y].push_back(make_pair(a[i].id, j));
			continue;
		}
		if (gety(y0) > gety(y1)) {
			swap(x0, x1);
			swap(y0, y1);
		}
		double k, b;
		k = (x1 - x0) / (y1 - y0);
		b = x1 - k * y1;  // x=ky+b
		if (k > 0) {
			int x = getx(x0), xm = getx(x1);
			int y = gety(y0), ym = gety(y1);
			for (; x >= xm && y <= ym; x--)
				for (; y <= ym; y++) {
					getleft_up(x, y, x0, y0);
					getright_down(x, y, x1, y1);
					if ((x0 - k * y0 - b) * (x1 - k * y1 - b) > 0) {
						y--;
						break;
					}
					else
						matrix[x][y].push_back(make_pair(a[i].id, j));
				}
		}
		else {
			int x = getx(x0), xm = getx(x1);
			int y = gety(y0), ym = gety(y1);
			for (; x <= xm && y <= ym; x++)
				for (; y <= ym; y++) {
					getleft_down(x, y, x0, y0);
					getright_up(x, y, x1, y1);
					if ((x0 - k * y0 - b) * (x1 - k * y1 - b) > 0) {
						y--;
						break;
					}
					else
						matrix[x][y].push_back(make_pair(a[i].id, j));
				}
		}
	}
}

inline double distance(double x, double y, pair<int, int> roadseg, double* w, double* j) {
	double x0 = a[roadseg.first].p[roadseg.second - 1].x, y0 = a[roadseg.first].p[roadseg.second - 1].y;
	double x1 = a[roadseg.first].p[roadseg.second].x, y1 = a[roadseg.first].p[roadseg.second].y;
	double cross = (x1 - x0) * (x - x0) + (y1 - y0) * (y - y0);
	if (cross <= 0) {
		*w = x0;
		*j = y0;
		return sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
	}
	double d = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
	if (cross >= d) {
		*w = x1;
		*j = y1;
		return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
	}
	double r = cross / d;
	*w = x0 + (x1 - x0) * r;
	*j = y0 + (y1 - y0) * r;
	return sqrt((x - *w) * (x - *w) + (y - *j) * (y - *j));
}

inline void FindMS(vector<vector<condidate>>& TMP, vector<pointT>& T) {
	vector<int>tmp;
	vector<PRE>pre(TMP.size());
	for (int j = 0; j < TMP[0].size(); j++)
		pre[0].F[j] = TMP[0][j].B;

	for (int i = 1; i < TMP.size(); i++) {
		double dis_o, cos, v;
		dis_o = sqrt((T[i].x - T[i - 1].x) * (T[i].x - T[i - 1].x) + (T[i].y - T[i - 1].y) * (T[i].y - T[i - 1].y));

		int Bk = 0;
		double fiB = 0, seB = 0;

		int tag = 0;
		double chongzhi;
		for (int j = 0; j < TMP[i].size(); j++) {
			double  F_B = 0, beita, fxs, median;
			int qian = 0;
			vector<double>dis_r(TMP[i - 1].size());
			vector<double>dis(TMP[i - 1].size());
			vector<double>F(TMP[i - 1].size());
			for (int k = 0; k < TMP[i - 1].size(); k++) {
				(PDJC(TMP[i - 1][k], TMP[i][j], F[k], dis_r[k],T[i].t-T[i-1].t));
				dis[k] = abs(dis_o - dis_r[k]);
			}
			for (int k = 0; k < TMP[i - 1].size(); k++) {
				F[k] *= Fxs * exp(-dis[k] * Fxs);
				F[k] *= pre[i - 1].F[k];
				if (F[k] > F_B) {
					F_B = F[k];
					qian = k;
				}
			}
			if (tag == 0 && j == 0 && F_B > 1e5) {
				tag = 1;
				chongzhi = F_B;
				F_B = 1;
			}
			else if (tag == 1) {
				F_B /= chongzhi;
			}
			F_B *= TMP[i][j].B;
			pre[i].F[j] = F_B;
			pre[i].qian[j] = qian;
		}
	}
	int qian = 0, F_B = 0, now = 0;
	for (int k = 0; k < SIZE; k++) {
		if (pre[pre.size() - 1].F[k] > F_B) {
			F_B = pre[pre.size() - 1].F[k];
			qian = pre[pre.size() - 1].qian[k];
			now = k;
		}
	}
	tmp.push_back(TMP[TMP.size() - 1][now].idfi);
	for (int i = TMP.size() - 2; i >= 0; i--) {
		tmp.push_back(TMP[i][qian].idfi);
		qian = pre[i].qian[qian];
	}
	reverse(tmp.begin(), tmp.end());
	ANS.push_back(tmp);
}

inline void GetCandidates(vector<pointT>& T) {
	vector<vector<condidate>> TMP;
	int x, y;
	double W, J, x1, x2, y1, y2, x0, y0;
	for (int k = 0; k < T.size(); k++) {
		vector<condidate> tmp, tmp2;
		W = T[k].x;
		J = T[k].y;
		x = getx(W);
		y = gety(J);
		for (int i = x - 1; i <= x + 1; i++)
			for (int j = y - 1; j <= y + 1; j++) {
				for (int c = 0; c < matrix[i][j].size(); c++) {
					if (distance(W, J, matrix[i][j][c], &x0, &y0) <= fiftymeter) {
						double dis = 10000 * ONEW * ((W - x0) * (W - x0) + (J - y0) * (J - y0)); //  "1000"替换"ONEW" 调节参数
						double b = 20 * Bxs * exp(-0.5 * dis / (sigma * sigma));   // "20"调节参数
						tmp.push_back(condidate(matrix[i][j][c], x0, y0, b));
					}
				}
			}
		sort(tmp.begin(), tmp.end(), cmp);
		for (int i = 1; i < tmp.size();) {
			if (tmp[i].idfi == tmp[i - 1].idfi && tmp[i].idse == tmp[i - 1].idse) {
				tmp.erase(tmp.begin() + i);
			}
			else i++;
		}
		for (int i = 0; i < SIZE && i < tmp.size(); i++)
			tmp2.push_back(tmp[i]);
		
		TMP.push_back(tmp2);
	}
	FindMS(TMP, T);
}


int main() {
	cin >> n;
	double x, y;
	for (int i = 0; i < n; i++) {
		cin >> a[i].id >> a[i].idbegin >> a[i].idend >> a[i].way >> a[i].wtype >> a[i].c;
		grid[a[i].idbegin].push_back(idbegend(a[i].id, 1));
		grid[a[i].idend].push_back(idbegend(a[i].id, 2));
		for (int j = 1; j <= a[i].c; j++) {
			cin >> x >> y;
			a[i].p.push_back(point(x, y));
		}
		a[i].Dx = a[i].p[a[i].p.size() - 1].x - a[i].p[0].x;
		a[i].Dy = a[i].p[a[i].p.size() - 1].y - a[i].p[0].y;
		a[i].len = sqrt(pow(a[i].p[a[i].p.size() - 1].x - a[i].p[0].x, 2) + pow(a[i].p[a[i].p.size() - 1].y - a[i].p[0].y, 2));
		addroad(i);
	}
	cin >> m;
	int t, ddl = 0;
	for (int i = 1; i <= m; i++) {
		vector<pointT>T;
		cin >> t;
		while (t != ddl) {
			cin >> x >> y;
			T.push_back(pointT(t, x, y));
			cin >> t;
		}
		GetCandidates(T);
		ddl++;
	}
	cout << ddl << "\n";
	for (int i = 0; i < ANS.size(); i++) {
		for (int j = 0; j < ANS[i].size(); j++)
			cout << ANS[i][j] << " ";
		cout << "\n";
	}
	return 0;
}
