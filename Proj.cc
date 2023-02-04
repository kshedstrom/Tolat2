#include <iostream>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include "Proj.h"

using namespace::std;

const double Proj::PI = 3.14159265358979;
const double Proj::RTOD = 57.2957795130823;
const double Proj::DTOR = .017453292519943;
const double Proj::DTRH = .008726646259972;
const double Proj::REarth = 6.3708e6;

static double
SIGN(double a1,double a2)
{
    return(a2 >= 0. ? fabs(a1) : -fabs(a1));
}

template <class T>
T min2(T a, T b)
{
    return ((a < b) ? a : b);
}

template <class T>
T max2(T a, T b)
{
    return ((a > b) ? a : b);
}

void
Proj::map_init()
{
    double tmp1, tmp2;
    switch (_proj) {
    case MERCATOR:
	if (fabs(_rota) <= 0.0001 && fabs(_lat) <= 0.0001) {
	    sino = 1.0;
	    coso = 0.0;
	    sinr = 0.0;
	    cosr = 1.0;
	    _proj = FAST_MERC;
	}
	else if (fabs(_rota) >= 179.9999 && fabs(_lat) <= 0.0001) {
	    sino = -1.0;
	    _lon = _lon + 180.0;
	    coso = 0.0;
	    sinr = 0.0;
	    cosr = 1.0;
	    _proj = FAST_MERC;
	} else {
	    tmp1 = _rota*DTOR;
	    tmp2 = _lat*DTOR;
	    sinr = sin(tmp1);
	    cosr = cos(tmp1);
	    sino = sin(tmp2);
	    coso = cos(tmp2);
	    double sint = coso*cosr;
	    double cost = sqrt(1. - sint*sint);
	    tmp1 = sinr / cost;
	    tmp2 = sino / cost;
	    _lon = _lon - atan2(tmp1, -cosr*tmp2) * RTOD;
	    sinr = tmp1 * coso;
	    cosr = - tmp2;
	    sino = sint;
	    coso = cost;
	}
	break;
    case CONIC:
	{
	double avg = 0.5*(_lat+_rota);
	sino = SIGN(1.0, avg);
	double chi1 = (90.0 - sino*_lat)*DTOR;
	if (_lat==_rota)
	    coso = cos(chi1);
	else {
	    double aa, bb, cc, dd;
	    double chi2 = (90.0 - sino*_rota)*DTOR;
	    aa = sin(chi1);
	    bb = sin(chi2);
	    cc = tan(0.5*chi1);
	    dd = tan(0.5*chi2);
	    coso = log(aa/bb)/log(cc/dd);
	}
	}
	break;
    case STEREOGRAPHIC:
	tmp1 = _rota*DTOR;
	tmp2 = _lat*DTOR;
	sinr = sin(tmp1);
	cosr = cos(tmp1);
	sino = sin(tmp2);
	coso = cos(tmp2);
	break;
    default:
    /* error, you didn't specify a projection we know */
	cerr << "ERROR in map_init -- no projection specified.\n";
	exit(1);
    }
}


void
Proj::map_tran(double ilat, double ilon, double& u, double& v)
{
    /* local variables */
    double tmp1, tmp2, cosa, sina, cosb, sinb;
    double chi, r, sinph, sinla, cosph, cosla, amin1;
    double tcos;

    /* diagnostics variables ********/
    double sign_arg1, sign_arg2;

  /*
   * set up u and v for the for the fast paths.
   *  u ==> longitude in degrees. -180.0 <= u <= 180.0
   *  v ==> latitude in degrees.
   */

    tmp1 = ilon - _lon;
    sign_arg1 = SIGN(180.0, tmp1+180.0);
    sign_arg2 = SIGN(180.0, 180.0-tmp1);
    u = tmp1 - sign_arg1 + sign_arg2;
    v = ilat;

    switch (_proj) {
  /* fast path cylindrical projections require that _lat==_rota==0.0 */
    case FAST_MERC:
	if (fabs(ilat) > 89.9999){
	    /* invisible or undefined projection */
	    perror("Undefined projection in map_translate; unable to do fast-path Mercator");
	    exit(1);
	}
	u = u*DTOR;
	v = log(tan((ilat+90.0)*DTRH));
	break;
    case CONIC:
	chi = 90.0 - sino*ilat;
	if (chi >= 179.9999){
	    perror("Undefined projection in map_translate; unable to do conformal conic");
	    exit(1);
	}
	r = pow(tan(DTRH*chi),coso);
	u = u * coso * DTOR;
	v = -r * sino * cos(u);
	u = r * sin(u);
	break;
    case STEREOGRAPHIC:
    case MERCATOR:
	tmp1 = u * DTOR;
	tmp2 = v * DTOR;
	sinph = sin(tmp1);
	sinla = sin(tmp2);
	cosph = cos(tmp1);
	cosla = cos(tmp2);
	tcos = cosla * cosph;

	amin1 = min2(1.0,  sinla*sino + tcos*coso);
	cosa  = max2(-1.0, amin1);
	sina = sqrt(1.0 - cosa*cosa);

	if (sina < 0.0001) {   /* error trap */
	    sina = 0.0;
	    if (cosa < 0.0){
		/*invisible or undefined projection */
		perror("Error in map_translate; cannot initialize stereographic projection");
		exit(1);
	    }
	    u = 0.0;
	    v = 0.0;
	    /* if this point is reached, we have bad initialization */
	    /* parameters bomb out and indicate this */
	    perror("Bad initial parameters in map_translate, no projection can be done");
	    exit(2);
	} /* end error trap */

	sinb = cosla*sinph/sina;
	cosb = (sinla*coso - tcos*sino)/sina;

	if (_proj == STEREOGRAPHIC) {
	    if (fabs(sina) < 0.0001)
		r = sina/2.0;
	    else
	        r = (1.0 - cosa)/sina;
	    u = r*(sinb*cosr + cosb*sinr);
	    v = r*(cosb*cosr - sinb*sinr);
	} /* end stereographic */
	else {
	    u = atan2(sinb*cosr + cosb*sinr, sinb*sinr - cosb*cosr);
	    v = log((1. + cosa)/sina);
	}
    } /* end switch */
} /* end function map_tran */

void
Proj::map_inv(double& ilat, double& ilon, double u, double v)
{
    /* local variables */
    double tmp1, tmp2, cosa, sina, cosb, sinb;
    double chi, r, sinph, sinla, cosph, cosla, amin1;
    double tcos;

    switch (_proj) {
  /* fast path cylindrical projections require that _lat==_rota==0.0 */
    case FAST_MERC:
//	if (fabs(u) > PI) {
//	    cerr << "inverse is not defined.\n";
//	    exit(2);
//	}
	ilat = 2. * RTOD * atan(exp(v)) - 90.;
	ilon = _lon + RTOD*u;
	break;
    case CONIC:
	r = sqrt(u*u + v*v);
	if (r < 1.e-6) {
	    ilat = sino*90.;
	    tmp1 = 0;
	    tmp2 = 1;
	} else {
	    ilat = sino * (90. - 2.*RTOD*atan(pow(r,(1./coso))));
	    tmp1 = u / r;
	    tmp2 = -sino * v / r;
	}
        ilon = _lon + RTOD*atan2(tmp1, tmp2) / coso;
//	if (fabs(ilon - _lon) > 180.) {
//	    cerr << "inverse is not defined.\n";
//	    exit(2);
//	}
	break;
    case STEREOGRAPHIC:
	r = sqrt(u*u + v*v);
	if (r < 1.e-6) {
	    sinb = 0.;
	    cosb = 1.;
	    sina = 0.;
	    cosa = 1.;
	} else {
	    sinb = (u*cosr - v*sinr) / r;
	    cosb = (u*sinr + v*cosr) / r;
	    sina = 2 * r / (1 + r*r);
	    cosa = (1. - r*r) / (1. + r*r);
	}
	sinph = sina * sinb;
	cosph = cosa * coso - sina * sino * cosb;
	cosla = sqrt(sinph*sinph + cosph*cosph);
	if (cosla != 0.) {
	    sinph = sinph / cosla;
	    cosph = cosph / cosla;
	}
	if (fabs(sino) > 1.e-6) {
	    sinla = (cosa - cosla * cosph * coso) / sino;
	} else {
	    sinla = sina * cosb;
	}
	ilat = RTOD * atan2(sinla, cosla);
	ilon = _lon + RTOD * atan2(sina*sinb, cosa*coso - sina*sino*cosb);
	break;
    case MERCATOR:
	{
	if (fabs(u) > PI) {
	    cerr << "inverse is not defined.\n";
	    exit(2);
	}
	sina = sin(PI - 2. * atan(exp(v)));
	cosa = cos(PI - 2. * atan(exp(v)));
	double sinu = sin(u);
	double cosu = cos(u);
	sinb = sinu*cosr + cosu*sinr;
	cosb = sinu*sinr - cosu*cosr;
	sinph = sina * sinb;
	cosph = cosa * coso - sina * sino * cosb;
	cosla = sqrt(sinph*sinph + cosph*cosph);
	if (cosla != 0.) {
	    sinph = sinph / cosla;
	    cosph = cosph / cosla;
	}
	if (fabs(sino) > 1.e-4) {
	    sinla = (cosa - cosla * cosph * coso) / sino;
	} else {
	    sinla = sina * cosb;
	}
	ilat = RTOD * atan2(sinla, cosla);
	ilon = _lon + RTOD * atan2(sina*sinb, cosa*coso - sina*sino*cosb);
	}
    } /* end switch */
    if (fabs(ilon) > 180.) ilon = ilon - SIGN(360., ilon);
} /* end function map_inv */
