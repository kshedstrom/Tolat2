#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <netcdf>
#include <math.h>
#include "Proj.h"

using namespace std;
using namespace netCDF;

int main(int argc, char** argv )
{

// Polar Stereographic - Arctic
//    double plon = 160.;
//    double rota = 0.;
//    double plat = 90.;
//    Proj p(Proj::STEREOGRAPHIC, plat, plon, rota);

// Polar Stereographic - Chukchi
//    double plon = -120.;
//    double rota = 0.;
//    double plat = 90.;
//    Proj p(Proj::STEREOGRAPHIC, plat, plon, rota);

// Cook Inlet
//    double plon = -95.;
//    double plon = -90.;
//    double plat = 55.;
//    double rota = 60.;
//    Proj p(Proj::CONIC, plat, plon, rota);

// Barrow canyon
//    double plat = 71.;
//    double rota = 73.;
//    double plon = -156.;
//    Proj p(Proj::CONIC, plat, plon, rota);

// Svalbard
//    double plon = 45.;
//    double plat = 78.;
//    double rota = 78.;
//    Proj p(Proj::CONIC, plat, plon, rota);

// CCS
//    double plon = -73.;
//    double rota = 40.;
//    double plat = 30.;
//    Proj p(Proj::CONIC, plat, plon, rota);
//
//    ???
//    double plon = -69.;
//    double rota = 45.;
//    double plat = 30.;
//    Proj p(Proj::CONIC, plat, plon, rota);

//  Beaufort
//    double plon = -162.;
//    double plon = -163.;
//    double rota = 70.;
//    double plat = 65.;
//    Proj p(Proj::CONIC, plat, plon, rota);

// NEP past
//    double plon = 180.;
//    double plon = -150.;
//    double plon = -139.;
//    double rota = 63.;
//    double plat = 57.;
// NEP7 + CGOA2.5
    double plon = -91.;
    double rota = 60.;
    double plat = 40.;
    Proj p(Proj::CONIC, plat, plon, rota);

// Caribbean
//    double plon = -80.;
// RCCS
//    double plon = -125.;
// Jinbo
//    double plon = 0.;
//    double rota = 0.;
//    double plat = 0.;
//    Proj p(Proj::MERCATOR, plat, plon, rota);

// Palau
//    double plon = 135.;
//    double rota = 0.;
//    double plat = 0.;
//    Proj p(Proj::MERCATOR, plat, plon, rota);

// Lombok Strait/SCS
//    double plon = 120.;
//    double rota = 0.;
//    double plat = 0.;
//    Proj p(Proj::MERCATOR, plat, plon, rota);

    double u1, u2, v1, v2;
    p.map_tran(plat+.5,plon,u2,v2);
    p.map_tran(plat-.5,plon,u1,v1);
    double udeg = sqrt((u2-u1)*(u2-u1) + (v2-v1)*(v2-v1));
    string gridfile;

    int npts;

    if (argc == 1) {
      cerr << "Name of netCDF file: ";
      cin >> gridfile;
    }
    else {
      gridfile = string(argv[1]);
    }

    NcFile nc(gridfile.c_str(), NcFile::write);
    NcDim d_xi_rho = nc.getDim("xi_rho");
    NcDim d_eta_rho = nc.getDim("eta_rho");
    NcDim d_xi_psi = nc.getDim("xi_psi");
    NcDim d_eta_psi = nc.getDim("eta_psi");

    int xi_rho = d_xi_rho.getSize();
    int eta_rho = d_eta_rho.getSize();
    int xi_psi = d_xi_psi.getSize();
    int eta_psi = d_eta_psi.getSize();

    NcVar x_psi = nc.getVar("x_psi");
    NcVar y_psi = nc.getVar("y_psi");
    NcVar lat_psi = nc.getVar("lat_psi");
    NcVar lon_psi = nc.getVar("lon_psi");
    double* xp = new double[xi_psi * eta_psi];
    double* yp = new double[xi_psi * eta_psi];
    double* lonp = new double[xi_psi * eta_psi];
    double* latp = new double[xi_psi * eta_psi];
    x_psi.getVar(xp);
    y_psi.getVar(yp);

    double u, v;
    for (int i=0; i < xi_psi*eta_psi; ++i)
        {
	  u = xp[i] * udeg * Proj::RTOD / Proj::REarth;
	  v = yp[i] * udeg * Proj::RTOD / Proj::REarth;
	  p.map_inv(latp[i], lonp[i], u, v);
        }

    vector<size_t> startp,countp;
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(eta_psi);
    countp.push_back(xi_psi);
    lat_psi.putVar(startp, countp, latp);
    lon_psi.putVar(startp, countp, lonp);

    NcVar x_rho = nc.getVar("x_rho");
    NcVar y_rho = nc.getVar("y_rho");
    NcVar lat_rho = nc.getVar("lat_rho");
    NcVar lon_rho = nc.getVar("lon_rho");
    double* xr = new double[xi_rho * eta_rho];
    double* yr = new double[xi_rho * eta_rho];
    double* lonr = new double[xi_rho * eta_rho];
    double* latr = new double[xi_rho * eta_rho];
    x_rho.getVar(xr);
    y_rho.getVar(yr);

     for (int i=0; i < xi_rho*eta_rho; ++i)
        {
	  u = xr[i] * udeg * Proj::RTOD / Proj::REarth;
	  v = yr[i] * udeg * Proj::RTOD / Proj::REarth;
	  p.map_inv(latr[i], lonr[i], u, v);
        }
    countp[0] = eta_rho;
    countp[1] = xi_rho;
    lat_rho.putVar(startp, countp, latr);
    lon_rho.putVar(startp, countp, lonr);

    NcVar x_u = nc.getVar("x_u");
    NcVar y_u = nc.getVar("y_u");
    NcVar lat_u = nc.getVar("lat_u");
    NcVar lon_u = nc.getVar("lon_u");
    double* xu = new double[xi_psi * eta_rho];
    double* yu = new double[xi_psi * eta_rho];
    double* lonu = new double[xi_psi * eta_rho];
    double* latu = new double[xi_psi * eta_rho];
    x_u.getVar(xu);
    y_u.getVar(yu);

    for (int i=0; i < xi_psi*eta_rho; ++i)
        {
  	  u = xu[i] * udeg * Proj::RTOD / Proj::REarth;
	  v = yu[i] * udeg * Proj::RTOD / Proj::REarth;
	  p.map_inv(latu[i], lonu[i], u, v);
        }
    countp[1] = xi_psi;
    lat_u.putVar(startp, countp, latu);
    lon_u.putVar(startp, countp, lonu);

    NcVar x_v = nc.getVar("x_v");
    NcVar y_v = nc.getVar("y_v");
    NcVar lat_v = nc.getVar("lat_v");
    NcVar lon_v = nc.getVar("lon_v");
    double* xv = new double[xi_rho * eta_rho];
    double* yv = new double[xi_rho * eta_rho];
    double* lonv = new double[xi_rho * eta_rho];
    double* latv = new double[xi_rho * eta_rho];
    x_v.getVar(xv);
    y_v.getVar(yv);

    for (int i=0; i < xi_rho*eta_rho; ++i)
    {
	u = xv[i] * udeg * Proj::RTOD / Proj::REarth;
	v = yv[i] * udeg * Proj::RTOD / Proj::REarth;
	p.map_inv(latv[i], lonv[i], u, v);
    }
    countp[0] = eta_psi;
    countp[1] = xi_rho;
    lat_v.putVar(startp, countp, latv);
    lon_v.putVar(startp, countp, lonv);

}
