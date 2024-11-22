//Author: Fukang Liu

#include <iostream>
#include <ctime>
#include <vector>
#include <thread>
#include <random>
using namespace std;


typedef unsigned long long u64;

struct TUPLE {
	int id;
	int x;
	int xp;
};

int rot8(int n, int s) {
	int numH = (n << s) & 0xff;
	int numL = (n >> (8 - s)) & 0xff;
	return (numH | numL);
}

int rot7(int n, int s) {
	int numH = (n << s) & 0x7f;
	int numL = (n >> (7 - s)) & 0x7f;
	return (numH | numL);
}

void TT(int arr[]) {
	for (int x = 0; x < 256; x++) {
		arr[x] = ((x + 1) * (x + 1) * (x + 1) - 1) % 257;
	}
}

void chi8(int arr[]) {
	for (int x = 0; x < 256; x++) {
		arr[x] = x ^ (rot8((~x)&0xff, 1) & rot8(x, 2) & rot8(x, 3));
		arr[x] = rot8(arr[x], 1);
	}
	//cout << "chi(255),chi(0):" << arr[255] << " " << arr[0] << endl;
	//cout << "chi(64),chi(63):" << arr[64] << " " << arr[63] << endl;
}

void chi7(int arr[]) {
	for (int x = 0; x < 128; x++) {
		arr[x] = x ^ (rot7((~x)&0x7f, 1) & rot7(x, 2));
		arr[x] = rot7(arr[x], 1);
	}
	//cout << "chi(127),chi(0):" << arr[127] << " " << arr[0] << endl;
}

void split31(u64 x, int arr[]) {
	for (int i = 0; i < 4; i++) {
		arr[i] = (x >> (i * 8)) & 0xff;
	}
}

void split64(u64 x, int arr[]) {
	for (int i = 0; i < 8; i++) {
		arr[i] = (x >> (i * 8)) & 0xff;
	}
}

u64 box31(u64 x,int arr8[],int arr7[]) {
	int arx[4], ary[4];
	u64 y = 0;
	for (int i = 0; i < 4; i++) {
		arx[i] = (x >> (i * 8)) & 0xff;
		if (i < 3)
			ary[i] = arr8[arx[i]];
		else
			ary[i] = arr7[arx[i]];
		y = y | ((u64(ary[i])) << (i * 8));
	}
	return y;
}

u64 box64(u64 x, int arr8[]) {
	int arx[8], ary[8];
	u64 y = 0;
	for (int i = 0; i < 8; i++) {
		arx[i] = (x >> (i * 8)) & 0xff;
		ary[i] = arr8[arx[i]];
		y = y | ((u64(ary[i])) << (i * 8));
	}
	return y;
}

u64 add64(u64 x, u64 y) {
	u64 modu = 18446744069414584321;//2^64 - 2^32 + 1
	u64 xh = (x >> 63) & 0x1;
	u64 xl = x & 0x7fffffffffffffff;
	u64 yh = (y >> 63) & 0x1;
	u64 yl = y & 0x7fffffffffffffff;
	u64 tmp = 0;
	if (xh + yh == 0) {
		return (xl + yl) % modu;
	}
	else if (xh + yh == 1) {
		if (xl + yl >= 0x8000000000000000) {
			tmp = (xl + yl) & 0x7fffffffffffffff;
			tmp = tmp + 0xffffffff;
			return tmp % modu;
		}
		else {
			return (0x8000000000000000 + xl + yl) % modu;
		}
	}
	else {
		return add64(0xffffffff, (xl + yl) % modu);
	}
}

u64 mul64(u64 x,u64 y) {
	u64 modu= 18446744069414584321;//2^64 - 2^32 + 1
	u64 xh = (x >> 32) & 0xffffffff;
	u64 xl = x & 0xffffffff;
	u64 yh = (y >> 32) & 0xffffffff;
	u64 yl = y & 0xffffffff;

	u64 tf1 = (xh * yh);
	u64 tf2 = (xh * yl);
	u64 tf3 = (xl * yh);
	u64 tf4 = (xl * yl);

	u64 s1 = add64(tf1, tf2);
	s1 = add64(s1, tf3);
	u64 sh = (s1 >> 32) & 0xffffffff;
	u64 sl = s1 & 0xffffffff;
	s1 = add64((0xffffffff * sh) % modu, (0x100000000 * sl) % modu);


	u64 s2 = modu - tf1;
	s2 = add64(s2, tf4);

	return add64(s1, s2);
}

void constructTAB(vector<vector<vector<int>>> &Tab,int size) {
	for (int carry = 0; carry < 2;  carry++) {
		for (int dx = 0; dx < size; dx++) {
			for (int x = 0; x < size; x++) {
				int xp = x + dx + carry;
				if (xp >= size) {
					Tab[carry][dx][x]= 1;
				}
				else {
					Tab[carry][dx][x] = 0;
				}
			}
		}
	}
}

void constructDDT(int j, int size, int arr[], vector<vector<vector<int>>>& Tab, vector<vector<vector<vector<TUPLE>>>>& DDT) {
	TUPLE tu;
	tu.id = 0;
	tu.x = 0;
	tu.xp = 0;
	int x, xp, y, yp,dx,dy;
	int cx1, cx2, cy1, cy2;
	for (int id = 0; id < 4; id++) {
		for (int t = 0; t < size; t++) {
			for (int tp = 0; tp < size; tp++) {
				//cout << id << " " << dx << " " << t << endl;
				cx1 = id / 2;
				cy1 = id % 2;
				x = t;
				xp = tp;
				if (j == 0 || j == 1) {
					//x = tp;
					//xp = t;
					dx = ((x + size * 2) - (cx1 + xp)) % size;
					cx2 = Tab[cx1][dx][xp];
				}
				else if (j == 2 || j == 3) {
					//x = t;
					//xp = tp;
					//x = tp;
					//xp = t;
					dx = ((xp + size * 2) - (cx1 + x)) % size;
					cx2 = Tab[cx1][dx][x];
				}
				tu.x = x;
				tu.xp = xp;
				y = arr[x];
				yp = arr[xp];
				if (j == 0 || j == 2) {
					dy = ((y + 2 * size) - (cy1 + yp)) % size;
					cy2 = Tab[cy1][dy][yp];
				}
				else if (j == 1 || j == 3) {
					dy = ((yp + 2 * size) - (cy1 + y)) % size;
					cy2 = Tab[cy1][dy][y];
				}
				tu.id = 2 * cx2 + cy2;
				DDT[id][dx][dy].push_back(tu);
			}
		}
	}
	/*for (int i = 0; i < 4; i++) {
		cout << "id="<<i<<" : " << DDT[i][0][0].size() << endl;
	}
	for (int i = 0; i < 4; i++) {
		cout << "id=" << i << " : ";
		for (int j = 0; j < DDT[i][0][0].size(); j++) {
			if (i!=0) {
				cout << DDT[i][0][0][j].id << " " << DDT[i][0][0][j].x << " " << DDT[i][0][0][j].xp << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
	*/
}

void constructTwoTables(int arr[], vector<vector<vector<int>>>& Tab, vector<vector<vector<vector<vector<TUPLE>>>>>& DDT, int size) {
	DDT.resize(4);
	for (int i = 0; i < 4; i++) {
		DDT[i].resize(4);
	}
	for (int tt = 0; tt < 4; tt++) {
		for (int i = 0; i < 4; i++) {
			DDT[tt][i].resize(size);
			for (int j = 0; j < size; j++) {
				DDT[tt][i][j].resize(size);
				for (int k = 0; k < size; k++) {
					DDT[tt][i][j][k].clear();
				}
			}
		}
	}

	Tab.resize(2);
	for (int i = 0; i < 2; i++) {
		Tab[i].resize(size);
		for (int j = 0; j < size; j++) {
			Tab[i][j].resize(size);
		}
	}

	constructTAB(Tab, size);
	constructDDT(0, size, arr, Tab, DDT[0]);
	constructDDT(1, size, arr, Tab, DDT[1]);
	constructDDT(2, size, arr, Tab, DDT[2]);
	constructDDT(3, size, arr, Tab, DDT[3]);
}




u64 getRandomNum(int len, u64 modu) {
	u64 x[8];
	u64 ex = 1 << 16;
	u64 val = 0;
	for (int i = 0; i < 4; i++) {
		x[i] = rand() % ex;
		val = val | (x[i] << (16 * i));
	}
	return val % modu;
}

u64 form31(u64 arr[]) {
	u64 val = 0;
	for (int i = 0; i < 4; i++) {
		val = val | (arr[i] << (8 * i));
	}
	return val;
}

u64 form64(u64 arr[]) {
	u64 val = 0;
	for (int i = 0; i < 8; i++) {
		val = val | (arr[i] << (8 * i));
	}
	return val;
}

void Re31(int j, int i, int id, 
	int h,
	int dw_arr[], int dz_arr[], 
	vector<u64>& wSol, 
	u64 w_arr[], u64 wp_arr[],
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT8,
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT7,
	u64 modu,int arr8[],int arr7[],u64 &checkTimes) {

	if (i == h) {
		u64 w = form31(w_arr);
		u64 wp = form31(wp_arr);
		if (w < modu && wp < modu && id==0) {
			wSol.push_back(w);
		}
	}
	else if(i<h) {
		int size = 0;
		int nextId = 0;
		if (i < h - 1) {//use DDT8
			size = DDT8[j][id][dw_arr[i]][dz_arr[i]].size();
			for (int j3 = 0; j3 < size; j3++) {
				nextId = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].id;
				w_arr[i] = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].x;
				wp_arr[i] = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].xp;
				checkTimes++;
				Re31(j, i + 1, nextId, h,  dw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, DDT7, modu, arr8, arr7,checkTimes);
				//cout << i << ", " << j << ": " << w_arr[i] << " " << z_arr[i] << "  size:" << size << endl;
			}
		}
		else {//use DDT7
			size = DDT7[j][id][dw_arr[i]][dz_arr[i]].size();
			for (int j3 = 0; j3 < size; j3++) {
				nextId = DDT7[j][id][dw_arr[i]][dz_arr[i]][j3].id;
				w_arr[i] = DDT7[j][id][dw_arr[i]][dz_arr[i]][j3].x;
				wp_arr[i] = DDT7[j][id][dw_arr[i]][dz_arr[i]][j3].xp;
				checkTimes++;
				Re31(j, i + 1, nextId, h,dw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, DDT7, modu, arr8, arr7,checkTimes);
			}
		}
	}

}

void Re64(int j, int i, int id,
	int h,
	int dw_arr[], int dz_arr[],
	vector<u64>& wSol,
	u64 w_arr[], u64 wp_arr[],
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT8,
	u64 modu, int arr8[],u64 &checkTimes) {

	if (i == h) {
		u64 w = form64(w_arr);
		u64 wp = form64(wp_arr);
		if (w < modu && wp < modu && id == 0) {
			wSol.push_back(w);
			wSol.push_back(wp);
		}
	}
	else if (i < h) {
		int size = 0;
		int nextId = 0;
		size = DDT8[j][id][dw_arr[i]][dz_arr[i]].size();
		for (int j3 = 0; j3 < size; j3++) {
			nextId = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].id;
			w_arr[i] = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].x;
			wp_arr[i] = DDT8[j][id][dw_arr[i]][dz_arr[i]][j3].xp;
			checkTimes++;
			Re64(j, i + 1, nextId, h, dw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, modu, arr8,checkTimes);
			//cout << i << ", " << j << ": " << w_arr[i] << " " << z_arr[i] << "  size:" << size << endl;
		}
	}

}

u64 findW_Monolith31(u64 dw, u64 dz, vector<u64>& wSol,u64 modu,
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT8,
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT7,
	int arr8[],int arr7[]) {
	u64 hdw = modu - dw;
	u64 hdz = modu - dz;
	int hdw_arr[4], hdz_arr[4], dw_arr[4], dz_arr[4];
	split31(hdw, hdw_arr);
	split31(hdz, hdz_arr);
	split31(dw, dw_arr);
	split31(dz, dz_arr);
	int h = 4;
	u64 w_arr[4], wp_arr[4];
	u64 checkTimes = 0;
	for (int j = 0; j < 4; j++) {
		int id = 0;
		int i = 0;
		if (j == 0) {
			Re31(j, i, id, h, hdw_arr, hdz_arr, wSol,w_arr,wp_arr,DDT8,DDT7,modu,arr8,arr7,checkTimes);
		}
		else if (j == 1) {
			Re31(j, i, id, h, hdw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, DDT7, modu, arr8, arr7, checkTimes);
		}
		else if (j == 2) {
			Re31(j, i, id, h, dw_arr, hdz_arr, wSol, w_arr, wp_arr, DDT8, DDT7, modu, arr8, arr7, checkTimes);
		}
		else if (j == 3) {
			Re31(j, i, id, h, dw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, DDT7, modu, arr8, arr7, checkTimes);
		}
	}
	return checkTimes;
}


u64 findW_64Version(u64 dw, u64 dz, vector<u64>& wSol, u64 modu,
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT8,
	int arr8[]) {
	u64 hdw = modu - dw;
	u64 hdz = modu - dz;
	int hdw_arr[8], hdz_arr[8], dw_arr[8], dz_arr[8];
	split64(hdw, hdw_arr);
	split64(hdz, hdz_arr);
	split64(dw, dw_arr);
	split64(dz, dz_arr);
	int h = 8;
	u64 w_arr[8], wp_arr[8];
	u64 checkTimes = 0;
	for (int j = 0; j < 4; j++) {
		int id = 0;
		int i = 0;
		if (j == 0) {
			Re64(j, i, id, h, hdw_arr, hdz_arr, wSol, w_arr, wp_arr, DDT8, modu, arr8,checkTimes);
		}
		else if (j == 1) {
			Re64(j, i, id, h, hdw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, modu, arr8,checkTimes);
		}
		else if (j == 2) {
			Re64(j, i, id, h, dw_arr, hdz_arr, wSol, w_arr, wp_arr, DDT8,  modu, arr8,checkTimes);
		}
		else if (j == 3) {
			Re64(j, i, id, h, dw_arr, dz_arr, wSol, w_arr, wp_arr, DDT8, modu, arr8,checkTimes);
		}
	}
	return checkTimes;
}


void getTabInfo(vector<vector<vector<vector<vector<TUPLE>>>>>& DDT, int size) {
	int nonzero[4] = { 0 };
	int carryzero[4] = { 0 };
	for (int t = 0; t < 4; t++) {
		for (int id = 0; id < 4; id++) {
			nonzero[id] = 0;
			carryzero[id] = 0;
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (DDT[t][id][i][j].size() > 0) {
						nonzero[id]++;
					}
					for (int z = 0; z < DDT[t][id][i][j].size(); z++) {
						if (DDT[t][id][i][j][z].id == 0) {
							carryzero[id]++;
						}
					}
				}
			}
			//cout << "j:"<<t<< " #non-empty entries:"<< nonzero[id] <<" #nextId=0:"<< carryzero[id] << endl;
		}
	}
}

void testMonolith31(u64 testTimes) {
	vector<vector<vector<int>>>Tab8; 
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT8;

	vector<vector<vector<int>>>Tab7; 
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT7;
	int arr8[256], arr7[128];
	chi8(arr8);
	chi7(arr7);

	constructTwoTables(arr8,Tab8, DDT8, 256);
	constructTwoTables(arr7,Tab7, DDT7, 128);
	getTabInfo(DDT8,256);
	getTabInfo(DDT7, 128);
	//cout << "start testing" << endl;
	srand(time(NULL));
	u64 modu = 2147483647;//2^31-1
	u64 success = 0;
	u64 totalSol = 0;
	vector<u64> wSol;
	u64 totalCheckTimes = 0;

	std::mt19937_64 mt;
	std::random_device rnd;
	mt.seed(rnd());

	for (int i = 0; i < testTimes; i++) {

		u64 dw = mt() % modu;
		u64 dz = mt() % modu;
		u64 w, wp, z, zp;
		

		wSol.clear();
		totalCheckTimes += findW_Monolith31(dw, dz, wSol, modu, DDT8, DDT7, arr8, arr7);

		for (int i = 0; i < wSol.size(); i++) {
			w = wSol[i];
			wp = (w + dw) % modu;
			z = box31(w, arr8, arr7);
			zp = box31(wp, arr8, arr7);
			if (dz != (zp + modu - z) % modu) {
				cout << "wrong catch:"<< w<<" "<<wp<<" "<<z << " " << zp << endl;
			}
		}
		if (wSol.size() > 0) {
			totalSol += wSol.size();
			success++;
		}
	}
	cout << "#total tests:" << dec << testTimes << endl;
	cout << "#valid input-output difference pairs:" << dec << success << endl;
	cout << "#total right pairs:" << dec << totalSol << endl;
	cout << "#average right pairs for valid input-output differences:" << double(totalSol) / success << endl;
	cout << "#average table lookups:" << totalCheckTimes / testTimes << endl;
}



void testMonolith64(u64 testTimes) {
	vector<vector<vector<int>>>Tab8;
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT8;

	int arr8[256];
	chi8(arr8);

	constructTwoTables(arr8, Tab8, DDT8, 256);
	getTabInfo(DDT8, 256);
	srand(time(NULL));
	u64 modu = 18446744069414584321;//2^64 - 2^32 + 1
	u64 success = 0;
	u64 totalSol = 0;
	vector<u64> wSol;
	u64 totalCheckTimes = 0;

	std::mt19937_64 mt;
	std::random_device rnd;
	mt.seed(rnd());

	for (int i = 0; i < testTimes; i++) {
		u64 dw = mt() % modu;
		u64 dz = mt() % modu;
		u64 w, wp, z, zp;

		wSol.clear();
		totalCheckTimes += findW_64Version(dw, dz, wSol, modu, DDT8, arr8);

		for (int i = 0; i < wSol.size()/2; i++) {
			w = wSol[2*i];
			wp = wSol[2 * i + 1];
			z = box64(w, arr8);
			zp = box64(wp, arr8);
			if (zp >= z && zp-z!=dz) {
				cout << "wrong catch:" << w << " " << wp << " " << z << " " << zp << endl;
			}
			else if (zp<z && modu-dz != z-zp) {
				cout << "wrong catch:" << w << " " << wp << " " << z << " " << zp << endl;
			}
		}
		if (wSol.size() > 0) {
			totalSol += wSol.size();
			success++;
		}
	}
	cout << "#total tests:" << dec << testTimes << endl;
	cout << "#valid input-output difference pairs:" << dec << success << endl;
	cout << "#total right pairs:" << dec << totalSol << endl;
	cout << "#average right pairs for valid input-output differences:" << double(totalSol) / success << endl;
	cout << "#average table lookups:" << totalCheckTimes / testTimes << endl;
}


void testTip5(u64 testTimes) {
	vector<vector<vector<int>>>Tab8;
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT8;

	int arr8[256];
	TT(arr8);

	constructTwoTables(arr8, Tab8, DDT8, 256);
	getTabInfo(DDT8, 256);
	//cout << "start testing" << endl;
	srand(time(NULL));
	u64 modu = 18446744069414584321;//2^64 - 2^32 + 1
	u64 success = 0;
	u64 totalSol = 0;
	vector<u64> wSol;
	u64 totalCheckTimes = 0;


	std::mt19937_64 mt;
	std::random_device rnd;
	mt.seed(rnd());

	for (int i = 0; i < testTimes; i++) {
		u64 dw = mt() % modu;
		u64 dz = mt() % modu;
		u64 w, wp, z, zp;

		wSol.clear();
		totalCheckTimes += findW_64Version(dw, dz, wSol, modu, DDT8, arr8);

		for (int i = 0; i < wSol.size() / 2; i++) {
			w = wSol[2 * i];
			wp = wSol[2 * i + 1];
			z = box64(w, arr8);
			zp = box64(wp, arr8);
			if (zp >= z && zp - z != dz) {
				cout << "wrong catch:" << w << " " << wp << " " << z << " " << zp << endl;
			}
			else if (zp < z && modu - dz != z - zp) {
				cout << "wrong catch:" << w << " " << wp << " " << z << " " << zp << endl;
			}
		}
		if (wSol.size() > 0) {
			totalSol += wSol.size();
			success++;
		}
	}
	cout << "#total tests:" << dec << testTimes << endl;
	cout <<"#valid input-output difference pairs:"<< dec << success<< endl;
	cout << "#total right pairs:"<< dec << totalSol << endl;
	cout << "#average right pairs for valid input-output differences:" << double(totalSol) / success << endl;
	cout << "#average table lookups:" << totalCheckTimes / testTimes << endl;
}



void solveDiffX(u64 dx[], u64 mx[], u64 modu, int boxNum) {
	for (int i = 0; i < boxNum; i++) {
		dx[i + 1] = (mx[i] * dx[0]) % modu;
	}
}

void solveDiffY(u64 dy[], u64 my[], u64 modu, int boxNum) {
	for (int i = 0; i < boxNum; i++) {
		u64 t1 = (my[2 * i] * dy[0]) % modu;
		u64 t2 = (my[2 * i + 1] * dy[boxNum + 1]) % modu;
		dy[i + 1] = (t1 + t2) % modu;
	}
}

void recursiveCheck(int threadNum, bool& found, int total, int layer, 
	u64 xv[], u64 xx[], u64 xxp[], 
	u64 dx[], u64 dy[], u64 times[], 
	vector<vector<vector<vector<vector<TUPLE>>>>> &DDT8,
	vector<vector<vector<vector<vector<TUPLE>>>>>& DDT7,
	int arr8[],int arr7[],
	u64 modu, int len) {

	vector<u64> pool;
	pool.clear();

	int i = layer;
	times[i]++;
	u64 qxx = ((xxp[i] * xxp[i]) % modu) + modu - ((xx[i] * xx[i]) % modu);
	qxx = qxx % modu;
	u64 dxx = (dy[i + 1] + modu - qxx) % modu;
	findW_Monolith31(dx[i + 1], dxx, pool, modu,DDT8, DDT7, arr8,arr7);//2nd S-box

	if (layer + 1 == total) {
		for (int j = 0; j < pool.size(); j++) {
			u64 x = pool[j];
			u64 xp = (x + dx[i + 1]) % modu;
			xx[i + 1] = box31(x,arr8,arr7);
			xxp[i + 1] = box31(xp,arr8,arr7);
			qxx = ((xxp[i + 1] * xxp[i + 1]) % modu) + modu - ((xx[i + 1] * xx[i + 1]) % modu);
			qxx = qxx % modu;
			xv[i + 1] = x;
			if (qxx == dy[i + 2]) {
				found = true;
				cout << dec<<"threadNum:"<<threadNum<<": found!" << endl;
				cout << "dx: ";
				for (int t = 0; t < total + 2; t++) {
					cout <<hex<< dx[t] << ", ";
				}
				cout << endl;
				cout << "x: ";
				for (int t = 0; t < total + 1; t++) {
					cout <<hex<< xv[t] << ", ";
				}
				cout << endl;
				cout << "S(x): ";
				for (int t = 0; box31(xv[t],arr8,arr7)==xx[t]&& t < total + 1; t++) {
					cout <<hex<< xx[t] << ", ";
				}
				cout << endl;
				cout << "S(x+dx): ";
				for (int t = 0; t < total + 1; t++) {
					cout <<hex<< xxp[t] << ", ";
				}
				cout << endl;
				cout << "dy: ";
				for (int t = 0; t < total + 3; t++) {
					cout <<hex<< dy[t] << ", ";
				}
				cout << endl;
				cout << "times:";
				for (int t = 0; t < total + 1; t++) {
					cout <<hex<< times[t] << " ";
				}
				cout << endl << endl;
				break;
			}
			times[i + 1]++;
		}
		pool.clear();
		return;
	}

	for (int j = 0; j < pool.size(); j++) {
		u64 x = pool[j];
		u64 xp = (x + dx[i + 1]) % modu;
		xx[i + 1] = box31(x, arr8, arr7);
		xxp[i + 1] = box31(xp,arr8,arr7);
		xv[i + 1] = x;
		recursiveCheck(threadNum, found, total, layer + 1, xv, xx, xxp, dx, dy, times, DDT8, DDT7, arr8,arr7, modu, len);
	}
	pool.clear();
}

void findDiffMonolith31(int threadNum) {
	const int boxNum = 8;
	u64 dx[boxNum + 1];
	u64 dxx[boxNum + 1];
	u64 dy[boxNum + 2];

	int len = 4;
	u64 modu = 2147483647;
	//srand(time(NULL)+threadNum);

	std::mt19937_64 mt;
	std::random_device rnd;
	mt.seed(rnd() + threadNum);

	//the linear eq systems are obtained using a sage script (sys.sage)
	u64 mxRes[8] = { 933109610, 1272471931, 1943874347,  913546822,  718398565, 1015931969, 1469263215, 1836315033 };
	//linear eq system for x after Gaussian elimination

	u64 myRes[8][2] = { 1861297064, 765893364,
	255954954, 1522741964,
	205871370,  146296496,
	833069128, 1062429990,
	1622731249, 1674883117,
	963374190, 1553944205,
	568839686,  749463949,
	1437408028, 1002402912 };
	//linear equation system for y after Gaussian elimination

	


	u64 mx[boxNum];
	u64 my[boxNum * 2];
	for (int i = 0; i < boxNum; i++) {
		mx[i] = modu - mxRes[i];
		my[2 * i] = modu - myRes[i][1];
		my[2 * i + 1] = modu - myRes[i][0];
	}

	u64 x[boxNum], xp[boxNum], xx[boxNum], xxp[boxNum],xv[boxNum];
	bool found = false;

	u64 times[boxNum] = { 0 };


	vector<vector<vector<int>>>Tab8;
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT8;
	vector<vector<vector<int>>>Tab7;
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT7;
	int arr8[256], arr7[128];
	chi8(arr8);
	chi7(arr7);
	constructTwoTables(arr8, Tab8, DDT8, 256);
	constructTwoTables(arr7, Tab7, DDT7, 128);

	//cout << "start searching..." << endl;

	while (!found) {
		x[0] = mt() % modu;
		xp[0] = mt() % modu;
		dx[0] = (xp[0] + modu - x[0]) % modu;

		//verify the results in the paper
		/*
		x[0] = 0x69737ce4;
		xp[0] =0x4a9f18c5;
		dx[0] = (xp[0] + modu - x[0]) % modu;
		*/

		xx[0] = box31(x[0],arr8,arr7);
		xxp[0] = box31(xp[0], arr8, arr7);
		dxx[0] = (xxp[0] + modu - xx[0]) % modu;

		dy[0] = dxx[0];

		solveDiffX(dx, mx, modu, boxNum);//solve linear eq
		dy[boxNum + 1] = dx[boxNum];
		solveDiffY(dy, my, modu, boxNum);//solve linear eq

		xv[0] = x[0];
		recursiveCheck(threadNum, found, boxNum - 1, 0, xv,xx, xxp, dx, dy, times, DDT8, DDT7, arr8,arr7, modu, len);
	}
}

void verifyMono31() {
	u64 x[8] = { 1769176292, 474569811, 2112689211, 1266486537, 1481086973, 1130812394, 1590835993, 1234430009};
	u64 dx[9] = { 1630247904, 2107149975, 1802005106, 1378426179, 350811560, 815160035, 1958855613, 1335208794, 645582796};
	u64 xp[8];
	u64 r[9], rp[9];

	u64 modu = 2147483647;
	int arr8[256], arr7[128];
	chi8(arr8);
	chi7(arr7);


	cout << "x:";
	for (int i = 0; i < 8; i++) {
		cout << hex << "0x" << x[i] << ", ";
	}
	cout << endl;

	cout << "x':";
	for (int i = 0; i < 8; i++) {
		xp[i] = (x[i] + dx[i]) % modu;
		r[i] = box31(x[i], arr8, arr7);
		rp[i] = box31(xp[i], arr8, arr7);
		cout <<"0x" << xp[i] << ", ";
	}
	cout << endl;

	cout << "dx:"<<dec;
	for (int i = 0; i < 24; i++) {
		if (i < 8)
			cout << "0x" << dx[i] << ", ";
		else if (i == 23)
			cout << "0x" << dx[8];
		else
			cout << 0 << ", ";
	}
	cout << endl;

	cout << "out of sbox:"<<hex;
	for (int i = 0; i < 8; i++) {
		cout << hex << "0x" << r[i] << ", ";
	}
	cout << endl;

	cout << "outp of sbox:";
	for (int i = 0; i < 8; i++) {
		cout << hex << "0x" << rp[i] << ", ";
	}
	cout << endl;

	cout << "dr:";
	for (int i = 0; i < 24; i++) {
		if (i < 8)
			cout << "0x" << ((modu-r[i])+rp[i])%modu << ", ";
		else if (i == 23)
			cout << "0x" << dx[8];
		else
			cout << 0 << ", ";
	}
	cout << endl;
	

	r[8] = 1;
	rp[8] = 1;

	u64 dy[10];
	u64 dyRes[10] = { 634041782, 915772315, 2038072678, 724312696, 1712880358, 198008050, 1958559545, 1759016090, 2000721589, 645582796 };
	dy[9] = dx[8];

	dy[0] = (rp[0] + modu - r[0]) % modu;
	for (int i = 1; i < 9; i++) {
		u64 tmp1 = ((rp[i - 1] * rp[i - 1]) % modu + rp[i]) % modu;
		u64 tmp2 = ((r[i - 1] * r[i - 1]) % modu + r[i]) % modu;
		dy[i] = (tmp1 + modu - tmp2) % modu;
	}
	for (int i = 0; i < 10; i++) {
		cout <<dec<< dy[i]-dyRes[i] << " ";
	}
	cout << endl;


	cout << "dy:"<<dec;
	for (int i = 0; i < 24; i++) {
		if (i < 9)
			cout << "0x" << dy[i] << ", ";
		else if (i == 23)
			cout << "0x" << dy[9];
		else
			cout << 0 << ", ";
	}
	cout << endl;
}

void findDiffTip5(int threadNum) {
	u64 modu = 18446744069414584321;//2^64 - 2^32 + 1
	u64 times = 0;
	//srand(time(NULL) + threadNum);

	std::mt19937_64 mt;
	std::random_device rnd;
	mt.seed(rnd() + threadNum);


	vector<vector<vector<int>>>Tab8;
	vector<vector<vector<vector<vector<TUPLE>>>>> DDT8;
	int arr8[256];
	TT(arr8);
	constructTwoTables(arr8, Tab8, DDT8, 256);
	vector<u64> wSol;
	wSol.clear();
	//getTabInfo(DDT8, 256);

	//m[4] and mInv[4] are obtained used sys.sage
	u64 mInv[4] = { 11271452075168550964,17842331458134305255,  9118427056206055769, 11160408412621113223 };
	//linear eq sys after Gaussian elimination

	u64 m[4] = { 13476799261015796163, 14726433080041406471, 12569196116954954328, 10117416465323008938 };
	//linear eq sys after Gaussian elimination

	u64 xd[4] = { 0 }, yd[4] = { 0 };

	bool find = false;
	u64 x1 = 0, x11, x12, x13, x14;
	u64 dx = 0,dy=0;
	u64 x2, x21, x22, x23, x24;
	u64 xSol[5] = { 0 }, xpSol[5]={0}, ySol[5] = { 0 };
	while (!find) {
		times++;
		x1 = mt() % modu;
		dx = mt() % modu;
		x2 = add64(x1, dx);

		//check the results in the paper efficiently.
		/*
		x1 = 0x9b348004e33e78c5;
		x2 = 0x658096b7fa6a0246;
		dx = modu - x1;
		dx = add64(x2, dx);//dx = x2 - x1 mod modu
		*/

		//compute x14=x1^7
		x11 = mul64(x1, x1);//x1^2
		x12 = mul64(x11, x11);//x1^4
		x13 = mul64(x11, x12);//x1^6
		x14 = mul64(x13, x1);//x1^7

		//compute x24=x2^7
		x21 = mul64(x2, x2);//x2^2
		x22 = mul64(x21, x21);//x2^4
		x23 = mul64(x21, x22);//x2^6
		x24 = mul64(x23, x2);//x2^7

		//compute output diff of the T-box
		dy = add64(x24, modu - x14);

		xSol[4] = x1;
		xpSol[4] = x2;

		//solve linear eq systems
		for (int i = 0; i < 4; i++) {
			xd[i] = mul64(dx, mInv[i]);
			yd[i] = mul64(dy, m[i]);
		}

		bool search = true;
		for (int i = 0; i < 4; i++) {
			wSol.clear();
			findW_64Version(xd[i], yd[i], wSol, modu, DDT8, arr8);
			if (wSol.size() == 0) {
				search = false;
				break;
			}
			else {
				xSol[i] = wSol[0];
				xpSol[i] = wSol[1];
			}
		}
		wSol.clear();
		if (search == true) {
			find = true;
		}
	}
	cout <<dec<< "ThreadNum:" << threadNum << endl;
	cout << hex << "times:" << times << endl;
	cout << "x:" << hex<<xSol[0] << " " << xSol[1] << " " << xSol[2] << " " << xSol[3] << " " << xSol[4] << endl;
	cout << "xp:" << hex<<xpSol[0] << " " << xpSol[1] << " " << xpSol[2] << " " << xpSol[3] << " " << xpSol[4] << endl;
	cout << "dx:" << hex<<xd[0] << ", " << xd[1] << ", " << xd[2] << ", " << xd[3]<< ", "<<dx << endl;
	cout << "dy:" << hex<<yd[0] << ", " << yd[1] << ", " << yd[2] << ", " << yd[3] << ", " << dy << endl;
}


int main() {
	u64 testTimes = 1000000;
	


	cout << "Test Monolith31 S-box:" << endl;
	testMonolith31(testTimes);
	cout << endl;

	cout << endl<<"Test Monolith64 S-box:" << endl;
	testMonolith64(testTimes);
	cout << endl;

	cout << endl << "Test Tip5 S-box:" << endl;
	testTip5(testTimes);
	cout << endl;
	
	//cout << "finish" << endl;

	vector<thread> myThreads;
	int threadNum = 16;
	cout<<"1->search differences for Tip5;"<<endl;
	cout << "2->search differences for Monolith31;" << endl;
	int cmd;
	cout << "Please input the command (input 1 or 2): ";
	cin >> cmd;

	cout << "Please input the number of used threads:";
	cin >> threadNum;
	for (int j = 0; j < threadNum; j++) {
		if (cmd == 2) {
			myThreads.push_back(thread(findDiffMonolith31, j));
		}
		else if (cmd == 1) {
			myThreads.push_back(thread(findDiffTip5, j));
		}
		else {
			cout << "wrong command, exit";
			return 0;
		}
	}
	for (thread& single : myThreads) {
		single.join();
	}
	return 1;
}