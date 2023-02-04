class Proj {
public:
    enum proj {MERCATOR, CONIC, STEREOGRAPHIC, FAST_MERC};
    Proj(proj p, double lat, double lon, double rot):
	_proj(p), _lat(lat), _lon(lon), _rota(rot) { map_init(); }
    void map_tran(double lat, double lon, double& u, double& v);
    void map_inv(double& lat, double& lon, double u, double v);
    proj type() const { return _proj; }
    double lat() const { return _lat; }
    double lon() const { return _lon; }
    double rota() const { return _rota; }

/* define some trig constants **********
     DTOR = PI/180 (degrees to radians)
     DTRH = DTOR/2 or PI/360
     RTOD = 180/PI (radians to degrees)
*****************************************/
    static const double PI;
    static const double RTOD;
    static const double DTOR;
    static const double DTRH;
    static const double REarth;

private:
    proj _proj;			// type of projection
    double _lat;		// central latitude
    double _lon;		// central longitude
    double _rota;		// rotation
    double sino, coso, sinr, cosr;

    void map_init();
};
