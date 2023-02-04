#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include "Proj.h"

using namespace::std;

int main()
{
//    double plon = 180.;
    double plon = -80.;
    double rota = 0.;
    double plat = 0.;
    Proj p(Proj::MERCATOR, plat, plon, rota);

    bool backwards = false;
    double lat, lon, udeg;
    double u1, u2, v1, v2, u, v;

    p.map_tran(plat+.5, plon, u2, v2);
    p.map_tran(plat-.5, plon, u1, v1);
    udeg = sqrt((u2-u1)*(u2-u1) + (v2-v1)*(v2-v1));
    cerr << "udeg = " << udeg << endl;

    int n = 0;
    vector<double> x;
    vector<double> y;
    while (cin >> lat >> lon) {
	if (lat < 99.) {
            ++n;
	    p.map_tran(lat, lon, u, v);
	    u = u * Proj::DTOR * Proj::REarth / udeg;
	    v = v * Proj::DTOR * Proj::REarth / udeg;
	    x.push_back(u);
	    y.push_back(v);
	} else {
            if (backwards) {
	        reverse(x.begin(), x.end());
	        reverse(y.begin(), y.end());
            }
            cout << n << endl;
            for (int i = 0; i < n; ++i) {
	        cout << "  " << setprecision(14) << x[i];
	        cout << "  " << setprecision(14) << y[i] << endl;
            }
	    n = 0;
	    x.resize(0);
	    y.resize(0);
	}
    }
    return 0;
}
