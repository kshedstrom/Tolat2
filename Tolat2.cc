#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <netcdf.hh>
#include <math.h>
#include "Proj.h"

using namespace::std;

int main(int argc, char** argv )
{

// Polar Stereographic - Arctic
    double plon = 160.;
    double rota = 0.;
    double plat = 90.;
    Proj p(Proj::STEREOGRAPHIC, plat, plon, rota);

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

// NEP
//    double plon = 180.;
//    double plon = -150.;
//    double plon = -139.;
//    double plon = -91.;
//    double rota = 63.;
//    double rota = 60.;
//    double plat = 57.;
//    double plat = 40.;
//    Proj p(Proj::CONIC, plat, plon, rota);

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

    NcFile nc(gridfile.c_str(), NcFile::Write);
    NcDim* d_xi_rho = nc.get_dim("xi_rho");
    NcDim* d_eta_rho = nc.get_dim("eta_rho");
    NcDim* d_xi_psi = nc.get_dim("xi_psi");
    NcDim* d_eta_psi = nc.get_dim("eta_psi");

    int xi_rho = d_xi_rho->size();
    int eta_rho = d_eta_rho->size();
    int xi_psi = d_xi_psi->size();
    int eta_psi = d_eta_psi->size();

    NcVar* x_psi = nc.get_var("x_psi");
    NcVar* y_psi = nc.get_var("y_psi");
    NcVar* lat_psi = nc.get_var("lat_psi");
    NcVar* lon_psi = nc.get_var("lon_psi");
    double* xp = new double[xi_psi * eta_psi];
    double* yp = new double[xi_psi * eta_psi];
    double* lonp = new double[xi_psi * eta_psi];
    double* latp = new double[xi_psi * eta_psi];
    x_psi->get(xp, x_psi->edges());
    y_psi->get(yp, y_psi->edges());

    double u, v;
    for (int i=0; i < xi_psi*eta_psi; ++i)
    {
	u = xp[i] * udeg * Proj::RTOD / Proj::REarth;
	v = yp[i] * udeg * Proj::RTOD / Proj::REarth;
	p.map_inv(latp[i], lonp[i], u, v);
    }
    lat_psi->put(latp, lat_psi->edges());
    lon_psi->put(lonp, lon_psi->edges());
    delete xp, yp, lonp, latp;

    NcVar* x_rho = nc.get_var("x_rho");
    NcVar* y_rho = nc.get_var("y_rho");
    NcVar* lat_rho = nc.get_var("lat_rho");
    NcVar* lon_rho = nc.get_var("lon_rho");
    double* xr = new double[xi_rho * eta_rho];
    double* yr = new double[xi_rho * eta_rho];
    double* lonr = new double[xi_rho * eta_rho];
    double* latr = new double[xi_rho * eta_rho];
    x_rho->get(xr, x_rho->edges());
    y_rho->get(yr, x_rho->edges());

    for (int i=0; i < xi_rho*eta_rho; ++i)
    {
	u = xr[i] * udeg * Proj::RTOD / Proj::REarth;
	v = yr[i] * udeg * Proj::RTOD / Proj::REarth;
	p.map_inv(latr[i], lonr[i], u, v);
    }
    lat_rho->put(latr, lat_rho->edges());
    lon_rho->put(lonr, lon_rho->edges());
    delete xr, yr, lonr, latr;

    NcVar* x_u = nc.get_var("x_u");
    NcVar* y_u = nc.get_var("y_u");
    NcVar* lat_u = nc.get_var("lat_u");
    NcVar* lon_u = nc.get_var("lon_u");
    double* xu = new double[xi_psi * eta_rho];
    double* yu = new double[xi_psi * eta_rho];
    double* lonu = new double[xi_psi * eta_rho];
    double* latu = new double[xi_psi * eta_rho];
    x_u->get(xu, x_u->edges());
    y_u->get(yu, x_u->edges());

    for (int i=0; i < xi_psi*eta_rho; ++i)
    {
	u = xu[i] * udeg * Proj::RTOD / Proj::REarth;
	v = yu[i] * udeg * Proj::RTOD / Proj::REarth;
	p.map_inv(latu[i], lonu[i], u, v);
    }
    lat_u->put(latu, lat_u->edges());
    lon_u->put(lonu, lon_u->edges());
    delete xu, yu, lonu, latu;

    NcVar* x_v = nc.get_var("x_v");
    NcVar* y_v = nc.get_var("y_v");
    NcVar* lat_v = nc.get_var("lat_v");
    NcVar* lon_v = nc.get_var("lon_v");
    double* xv = new double[xi_rho * eta_rho];
    double* yv = new double[xi_rho * eta_rho];
    double* lonv = new double[xi_rho * eta_rho];
    double* latv = new double[xi_rho * eta_rho];
    x_v->get(xv, x_v->edges());
    y_v->get(yv, x_v->edges());

    for (int i=0; i < xi_rho*eta_rho; ++i)
    {
	u = xv[i] * udeg * Proj::RTOD / Proj::REarth;
	v = yv[i] * udeg * Proj::RTOD / Proj::REarth;
	p.map_inv(latv[i], lonv[i], u, v);
    }
    lat_v->put(latv, lat_v->edges());
    lon_v->put(lonv, lon_v->edges());
    delete xv, yv, lonv, latv;

}
