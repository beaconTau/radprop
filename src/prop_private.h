
/** These are just headers for the definitions in itwom3.0.cpp 
 *
 *  see that file for documentation 
 * */ 

// itwomv3 
void point_to_point(double elev[], double tht_m, double rht_m, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double frq_mhz, int radio_climate, int pol, double conf, double rel, double &dbloss, char *strmode, int &errnum,
    double clutter_refractivity = 1000, double clutter_height = 22.5, double clutter_density = 1,int mode_var = 1,  double delta_h_diff = 0); 

//itm /l-r 
void point_to_point_ITM(double elev[], double tht_m, double rht_m, double eps_dielect, double sgm_conductivity, double eno_ns_surfref, double frq_mhz, int radio_climate, int pol, double conf, double rel, double &dbloss, char *strmode, int &errnum); 
