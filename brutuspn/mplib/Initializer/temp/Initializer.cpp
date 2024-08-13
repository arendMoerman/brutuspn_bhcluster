#include "Initializer.h"

vector<mpreal> Initializer::generate(int N, string config, vector<string> par) {
  if(config == "file") {
    vector<mpreal> data;
    bool data_obtained = get_from_file(data, par[0]);
    if(data_obtained) return data;
  }
  if(config == "file_list") {  
    vector<mpreal> data;
    bool data_obtained = get_from_file_list(data, par[0], atoi(par[1].c_str()));
    if(data_obtained) return data;
  }
  if(config == "continue_from_file") {  
    vector<mpreal> data;
    bool data_obtained = continue_from_file(data, par[0]);
    if(data_obtained) return data;
  }

  if(N == 2) {
    if(config == "circle") return get_circle( atof(par[0].c_str()) );
    if(config == "binary_random") return get_binary();
    if(config == "binary") return get_binary(atof(par[0].c_str()), atof(par[1].c_str()));
    if(config == "parabola") return get_parabola(atof(par[0].c_str()));
    if(config == "hyperbola") return get_hyperbola(atof(par[0].c_str()), atof(par[1].c_str()));
    if(config == "radial") return get_radial( atof(par[0].c_str()), atof(par[1].c_str()) );
  }
  if(N == 3) {
    if(config == "figure8") return get_figure8();
    if(config == "figure8_precision") return get_figure8( atoi(par[0].c_str()) );
    if(config == "figure8_perturbed") return get_figure8_perturbed( atof(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()), atoi(par[3].c_str()) );

    if(config == "triangle") return get_triangle();
    if(config == "triangle_precision") return get_triangle( atoi(par[0].c_str()) );

    if(config == "butterfly") return get_butterfly();
    if(config == "butterfly_precision") return get_butterfly( atoi(par[0].c_str()) );

    if(config == "pythagorean") return get_pythagorean();
    if(config == "sitnikov") return get_sitnikov(atof(par[0].c_str()), atof(par[1].c_str()));
    if(config == "sitnikov_angle") return get_sitnikov(atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));

    if(config == "hierarchical_triple") return get_hierarchical_triple();

    if(config == "plummer_mass_ratio") return get_N3_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
    if(config == "cold_plummer_mass_ratio") return get_N3_cold_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
    if(config == "democratic_plummer_mass_ratio") return get_N3_democratic_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
    if(config == "democratic_cold_plummer_mass_ratio") return get_N3_democratic_cold_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
  }
  if(N == 4) {
    if(config == "plummer_mass_ratio") return get_N4_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
    if(config == "democratic_plummer_mass_ratio") return get_N4_democratic_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
    if(config == "democratic_cold_plummer_mass_ratio") return get_N4_democratic_cold_plummer_mass_ratio(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()));
  }
  if(N == 10) {
    if(config == "solar_system") return get_solar_system();
  }

  if(config == "full_solar_system") {
    return get_full_solar_system(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "sun_torus") {
    return get_sun_torus(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "sun_jupiter_torus") {
    return get_sun_jupiter_torus(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "sun_mars_jupiter_torus") {
    return get_sun_mars_jupiter_torus(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "sun_mars_jovian_torus") {
    return get_sun_mars_jovian_torus(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }

  if(config == "plummer") {
    return get_plummer(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "cold_plummer") {
    return get_cold_plummer(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }

  if(config == "democratic_plummer") {
    return get_democratic_plummer(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }
  if(config == "democratic_cold_plummer") {
    return get_democratic_cold_plummer(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()));
  }

  if(config == "bh_plummer") {
    return get_bh_plummer(N, atof(par[0].c_str()), atof(par[1].c_str()), atof(par[2].c_str()), atof(par[3].c_str()));
  }
  if(config == "bh_solar_systems") {
    int numPlanets = atoi(par[0].c_str());
    int numBH = 1;
    int numStar = (N-numBH)/(numPlanets+1);
    return get_bh_solar_systems(numBH+numStar, numPlanets, atof(par[1].c_str()), atof(par[2].c_str()), atoi(par[3].c_str()), atoi(par[4].c_str()), atoi(par[5].c_str()));
  }
  if(config == "plummer_solar_systems") {
    int numPlanets = atoi(par[0].c_str());
    int numStar = N/(1+numPlanets);
    return get_plummer_solar_systems(numStar, numPlanets, atof(par[1].c_str()), atoi(par[2].c_str()), atoi(par[3].c_str()), atoi(par[4].c_str()));
  }

  cerr << "Invalid initial condition!" << endl;
  exit(1);
}

// Tools
void Initializer::set_random_generator(int seed, int pivot) {
  random.set_seed(seed);
  random.set_pivot(pivot);
}

void Initializer::centralize(vector<mpreal> &data) {
  int N = data.size()/7;
  vector<mpreal> rcm = diagnostics.get_rcm(data);
  vector<mpreal> vcm = diagnostics.get_vcm(data);
  for(int i=0; i<N; i++) {
    data[i*7+1] -= rcm[0];
    data[i*7+2] -= rcm[1];
    data[i*7+3] -= rcm[2];
    data[i*7+4] -= vcm[0];
    data[i*7+5] -= vcm[1];
    data[i*7+6] -= vcm[2];
  }  
}
void Initializer::normalize_to_nbody_units(vector<mpreal> &data) {
  int N = data.size()/7;/*
  if(N>1) {
    mpreal M = diagnostics.get_mass(data);
    mpreal Cm = "1.0"/M;
    for(int i=0; i<N; i++) {
      data[i*7] *= Cm;
    }
    mpreal r_vir = diagnostics.get_virial_radius(data);
    mpreal sigma = diagnostics.get_velocity_disperion(data);
    mpreal Cr = "1.0"/r_vir;
    mpreal Cv = "1.0"/sqrt("2.0")/sigma;
    for(int i=0; i<N; i++) {
      data[i*7+1] *= Cr;
      data[i*7+2] *= Cr;
      data[i*7+3] *= Cr;
      data[i*7+4] *= Cv;
      data[i*7+5] *= Cv;
      data[i*7+6] *= Cv;
    }   
  }*/
}
void Initializer::normalize_to_boekholt_units(vector<mpreal> &data) {
  int N = data.size()/7;
/*  mpreal M = (mpreal)diagnostics.get_mass(data);
  mpreal Cm = "1.0"/M;
  for(int i=0; i<N; i++) {
    data[i*7] *= Cm;
  }
  mpreal r_har = (mpreal)diagnostics.get_harmonic_radius(data);
  mpreal Cr = "1.0"/r_har;
  for(int i=0; i<N; i++) {
    data[i*7+1] *= Cr;
    data[i*7+2] *= Cr;
    data[i*7+3] *= Cr;
  }  

  mpreal Cv = sqrt(Cm/Cr);
  for(int i=0; i<N; i++) {
    data[i*7+4] *= Cv;
    data[i*7+5] *= Cv;
    data[i*7+6] *= Cv;
  }  */
}
void Initializer::rescale_position(vector<mpreal> &data, mpreal Cr) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    data[i*7+1] *= Cr;
    data[i*7+2] *= Cr;
    data[i*7+3] *= Cr;
  } 
}
void Initializer::rescale_mass_position(vector<mpreal> &data, mpreal Cm, mpreal Cr) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    data[i*7] *= Cm;
  }
  for(int i=0; i<N; i++) {
    data[i*7+1] *= Cr;
    data[i*7+2] *= Cr;
    data[i*7+3] *= Cr;
  }  
  mpreal Cv = sqrt(Cm/Cr);
  for(int i=0; i<N; i++) {
    data[i*7+4] *= Cv;
    data[i*7+5] *= Cv;
    data[i*7+6] *= Cv;
  }  
}
void Initializer::make_cold(vector<mpreal> &data) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    data[i*7+4] = "0";
    data[i*7+5] = "0";
    data[i*7+6] = "0";
  } 
}
void Initializer::rotate_x(vector<mpreal> &data, mpreal angle) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal xn  = data[i*7+1];
    mpreal yn  = cos(angle)*data[i*7+2]-sin(angle)*data[i*7+3];
    mpreal zn  = sin(angle)*data[i*7+2]+cos(angle)*data[i*7+3];
    mpreal vxn = data[i*7+4];
    mpreal vyn = cos(angle)*data[i*7+5]-sin(angle)*data[i*7+6];
    mpreal vzn = sin(angle)*data[i*7+5]+cos(angle)*data[i*7+6];
    data[i*7+1] = xn;
    data[i*7+2] = yn;
    data[i*7+3] = zn;
    data[i*7+4] = vxn;
    data[i*7+5] = vyn;
    data[i*7+6] = vzn;
  }
}
void Initializer::rotate_y(vector<mpreal> &data, mpreal angle) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal xn  = cos(angle)*data[i*7+1]+sin(angle)*data[i*7+3];
    mpreal yn  = data[i*7+2];
    mpreal zn  = -sin(angle)*data[i*7+1]+cos(angle)*data[i*7+3];
    mpreal vxn = cos(angle)*data[i*7+4]+sin(angle)*data[i*7+6];
    mpreal vyn = data[i*7+5];
    mpreal vzn = -sin(angle)*data[i*7+4]+cos(angle)*data[i*7+6];
    data[i*7+1] = xn;
    data[i*7+2] = yn;
    data[i*7+3] = zn;
    data[i*7+4] = vxn;
    data[i*7+5] = vyn;
    data[i*7+6] = vzn;
  }
}
void Initializer::rotate_z(vector<mpreal> &data, mpreal angle) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    mpreal xn  = cos(angle)*data[i*7+1]-sin(angle)*data[i*7+2];
    mpreal yn  = sin(angle)*data[i*7+1]+cos(angle)*data[i*7+2];
    mpreal zn  = data[i*7+3];
    mpreal vxn = cos(angle)*data[i*7+4]-sin(angle)*data[i*7+5];
    mpreal vyn = sin(angle)*data[i*7+4]+cos(angle)*data[i*7+5];
    mpreal vzn = data[i*7+6];
    data[i*7+1] = xn;
    data[i*7+2] = yn;
    data[i*7+3] = zn;
    data[i*7+4] = vxn;
    data[i*7+5] = vyn;
    data[i*7+6] = vzn;
  }
}
void Initializer::flip_z(vector<mpreal> &data) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    data[i*7+3] *= "-1";
    data[i*7+6] *= "-1";
  } 
}
void Initializer::flip_v(vector<mpreal> &data) {
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    data[i*7+4] *= "-1";
    data[i*7+5] *= "-1";
    data[i*7+6] *= "-1";
  } 
}

mpreal Initializer::get_random_mass_ratio() {
  // f(m1/m2)=U(1,10)
  mpreal m_ratio = random.uniform(1.0, 10.0);
  return m_ratio;
}
mpreal Initializer::get_random_eccentricity() {
  // f(e)=2e
  bool accept = false;
  while(!accept) {
    mpreal trial = (mpreal)random.uniform(0.0, 1.0);
    mpreal p_accept = "2"*trial;
    mpreal p = "2"*(mpreal)random.uniform(0.0, 1.0);
    if(p < p_accept) {
      return trial;
    }
  } 
}

// From file
bool Initializer::get_from_file(vector<mpreal> &data, string file) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<double> > file_data = data_handler.get_data();
    vector<double> data0 = file_data[ 0 ];
    for(int i=0; i<data0.size(); i++) data.push_back( (mpreal)data0[i] ); 
    return true;
  }
}
bool Initializer::get_from_file_list(vector<mpreal> &data, string file, int index) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<double> > file_data = data_handler.get_data();
    vector<double> data0 = file_data[ index ];
    for(int i=0; i<data0.size(); i++) data.push_back( (mpreal)data0[i] ); 
    return true;
  }
}
bool Initializer::continue_from_file(vector<mpreal> &data, string file) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<double> > file_data = data_handler.get_data();
    vector<double> data0 = file_data[ file_data.size()-1 ];
    for(int i=0; i<data0.size(); i++) data.push_back( (mpreal)data0[i] ); 
    return true;
  }
}

// N=2
vector<mpreal> Initializer::get_circle(mpreal m_ratio) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1.0"/m_ratio;
  mpreal m2 = "1.0"/("1.0"+m_ratio);
  mpreal m1 = "1.0"-m2;

  mpreal a = "1";
  mpreal v = "1";

  data.push_back(m1);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back(a);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back(v);
  data.push_back("0");

  centralize(data);
}
vector<mpreal> Initializer::get_binary() {
  vector<mpreal> data;

  mpreal m_ratio = get_random_mass_ratio();
  mpreal eccentricity = get_random_eccentricity();  

  if(m_ratio < "1") m_ratio = "1.0"/m_ratio;
  mpreal m2 = "1.0"/("1.0"+m_ratio);
  mpreal m1 = "1.0"-m2;

  mpreal ra = "1"+eccentricity;
  mpreal rp = "1"-eccentricity;

  mpreal v_aph = "1";
  if(eccentricity != "0") {
    v_aph = sqrt( "2.0" * ("1.0"/rp-"1.0"/ra) / (ra*ra/(rp*rp)-"1.0") );
  }

  data.push_back(m1);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back(ra);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back(v_aph);
  data.push_back("0");

  centralize(data);

  return data;
}
vector<mpreal> Initializer::get_binary(mpreal m_ratio, mpreal eccentricity) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1.0"/m_ratio;
  mpreal m2 = "1.0"/("1.0"+m_ratio);
  mpreal m1 = "1.0"-m2;

  mpreal ra = "1"+eccentricity;
  mpreal rp = "1"-eccentricity;

  mpreal v_aph = "1";
  if(eccentricity != "0") {
    v_aph = sqrt( "2.0" * ("1.0"/rp-"1.0"/ra) / (ra*ra/(rp*rp)-"1.0") );
  }

  data.push_back(m1);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back(ra);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back(v_aph);
  data.push_back("0");

  centralize(data);

  return data;
}

vector<mpreal> Initializer::get_parabola(mpreal m_ratio) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1.0"/m_ratio;
  mpreal m2 = "1.0"/("1.0"+m_ratio);
  mpreal m1 = "1.0"-m2;

  mpreal rp = "2.0";
  mpreal vp = "1.0";

  data.push_back(m1);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back("0");
  data.push_back(-rp);
  data.push_back("0");
  data.push_back(vp);
  data.push_back("0");
  data.push_back("0");

  centralize(data);  

  return data;
}
vector<mpreal> Initializer::get_hyperbola(mpreal m_ratio, mpreal eccentricity) {
  vector<mpreal> data;

  if(m_ratio > "1") m_ratio = "1.0"/m_ratio;
  mpreal m1 = "1.0"/("1.0"+m_ratio);
  mpreal m2 = "1.0"-m1;

  if(eccentricity <= "1") eccentricity = 1.0 + 1e-12;
  mpreal e = "0.5";

  mpreal vp = sqrt( "2.0"*(e + "1.0"/("1.0"-eccentricity)) );
  mpreal rp = "2.0" / (vp*vp - "1");

  data.push_back(m1);
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back("0");
  data.push_back(-rp);
  data.push_back("0");
  data.push_back(vp);
  data.push_back("0");
  data.push_back("0");  

  centralize(data);

  return data;
}

vector<mpreal> Initializer::get_radial(mpreal m_ratio, mpreal energy) {
  vector<mpreal> data;

  if(m_ratio > "1") m_ratio = "1.0"/m_ratio;
  mpreal m1 = "1.0"/("1.0"+m_ratio);
  mpreal m2 = "1.0"-m1;

  mpreal x1 = "0";
  mpreal x2 = "1e-6";

  mpreal vx1 = "0";
  mpreal vx2 = sqrt("2.0"*(energy+"1.0"/x2));

  data.push_back(m1);
  data.push_back(x1);
  data.push_back("0");
  data.push_back("0");
  data.push_back(vx1);
  data.push_back("0");
  data.push_back("0");

  data.push_back(m2);
  data.push_back(x2);
  data.push_back("0");
  data.push_back("0");
  data.push_back(vx2);
  data.push_back("0");
  data.push_back("0");

  centralize(data);

  return data;
}

// N=3
vector<mpreal> Initializer::get_figure8() {
  vector<mpreal> data;
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-0.93240737");
  data.push_back("-0.86473146");
  data.push_back("0");

  data.push_back("1");
  data.push_back("0.9700436");
  data.push_back("-0.24308753");
  data.push_back("0");
  data.push_back("0.466203685");
  data.push_back("0.43236573");
  data.push_back("0");

  data.push_back("1");
  data.push_back("-0.9700436");
  data.push_back("0.24308753");
  data.push_back("0");
  data.push_back("0.466203685");
  data.push_back("0.43236573");
  data.push_back("0");
  return data;
}
vector<mpreal> Initializer::get_figure8(int numDigits) {
  vector<mpreal> data;

  if(numDigits < 1) numDigits = 1;
  if(numDigits > 9) numDigits = 9;

  mpreal precision = pow("10", -numDigits);

  string aa = "0.9324073700";
  string bb = "0.8647314600";
  string cc = "0.9700436000";
  string dd = "0.2430875300";
  string ee = "0.4662036850";
  string ff = "0.4323657300";

  string a1 = "0.";
  string b1 = "0.";
  string c1 = "0.";
  string d1 = "0.";
  string e1 = "0.";
  string f1 = "0.";
  for(int i=0; i<numDigits; i++) {
    a1 += aa[2+i];
    b1 += bb[2+i];
    c1 += cc[2+i];
    d1 += dd[2+i];
    e1 += ee[2+i];
    f1 += ff[2+i];
  }
  string a2 = "0.";
  string b2 = "0.";
  string c2 = "0.";
  string d2 = "0.";
  string e2 = "0.";
  string f2 = "0.";
  for(int i=0; i<numDigits+1; i++) {
    a2 += aa[2+i];
    b2 += bb[2+i];
    c2 += cc[2+i];
    d2 += dd[2+i];
    e2 += ee[2+i];
    f2 += ff[2+i];
  }

  mpreal a3 = a1;
  mpreal b3 = b1;
  mpreal c3 = c1;
  mpreal d3 = d1;
  mpreal e3 = e1;
  mpreal f3 = f1;
  mpreal a4 = a2;
  mpreal b4 = b2;
  mpreal c4 = c2;
  mpreal d4 = d2;
  mpreal e4 = e2;
  mpreal f4 = f2;

  mpreal diffa = a4-a3;
  mpreal diffb = b4-b3;
  mpreal diffc = c4-c3;
  mpreal diffd = d4-d3;
  mpreal diffe = e4-e3;
  mpreal difff = f4-f3;

  if(diffa > precision/"2"-precision/"20") a3 += precision;
  if(diffb > precision/"2"-precision/"20") b3 += precision;
  if(diffc > precision/"2"-precision/"20") c3 += precision;
  if(diffd > precision/"2"-precision/"20") d3 += precision;
  if(diffe > precision/"2"-precision/"20") e3 += precision;
  if(difff > precision/"2"-precision/"20") f3 += precision;

  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-1"*a3);
  data.push_back("-1"*b3);
  data.push_back("0");

  data.push_back("1");
  data.push_back(c3);
  data.push_back("-1"*d3);
  data.push_back("0");
  data.push_back(e3);
  data.push_back(f3);
  data.push_back("0");

  data.push_back("1");
  data.push_back("-1"*c3);
  data.push_back(d3);
  data.push_back("0");
  data.push_back(e3);
  data.push_back(f3);
  data.push_back("0");

  centralize(data);

  return data;
}
vector<mpreal> Initializer::get_figure8_perturbed(mpreal dr, int seed, int pivot, int index) {
  vector<mpreal> data;

  set_random_generator(seed, pivot);

  for(int i=0; i<=index; i++) {
    data = get_figure8();  
    mpreal p = (random.uniform(0.0, 4.0)-2.0)/(dr/10.0);
    data[3] += dr+p;
  }

  centralize(data);

  return data;
}

vector<mpreal> Initializer::get_triangle() {
  vector<mpreal> data;

  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-0.5");
  data.push_back(sqrt("0.75"));
  data.push_back("0");

  data.push_back("1");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-0.5");
  data.push_back("-1"*sqrt("0.75"));
  data.push_back("0");

  data.push_back("1");
  data.push_back("0.5");
  data.push_back(sqrt("0.75"));
  data.push_back("0");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");

  centralize(data);

  return data;
}
vector<mpreal> Initializer::get_triangle(int numDigits) {
  vector<mpreal> data;

  if(numDigits < 1) numDigits = 1;
  if(numDigits > 15) numDigits = 15;

  mpreal precision = pow("10", -numDigits);

  string cc = "0.8660254037844386";

  string c1 = "0.";
  for(int i=0; i<numDigits; i++) {
    c1 += cc[2+i];
  }
  string c2 = "0.";
  for(int i=0; i<numDigits+1; i++) {
    c2 += cc[2+i];
  }

  mpreal c3 = atof(c1.c_str());
  mpreal c4 = atof(c2.c_str());

  mpreal diffc = c4-c3;

  if(diffc > precision/"2"-precision/"20") c3 += precision;

  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-0.5");
  data.push_back(c3);
  data.push_back("0");

  data.push_back("1");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-0.5");
  data.push_back(-c3);
  data.push_back("0");

  data.push_back("1");
  data.push_back("0.5");
  data.push_back(c3);
  data.push_back("0");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");

  centralize(data);

  return data;
}

vector<mpreal> Initializer::get_butterfly() {
  vector<mpreal> data;

  mpreal a = "0.306892758965492";
  mpreal b = "0.125506782829762";

  data.push_back("1");
  data.push_back("-1");
  data.push_back("0");
  data.push_back("0");
  data.push_back(a);
  data.push_back(b);
  data.push_back("0");

  data.push_back("1");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back(a);
  data.push_back(b);
  data.push_back("0");

  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-2"*a);
  data.push_back("-2"*b);
  data.push_back("0");

  return data;
}
vector<mpreal> Initializer::get_butterfly(int numDigits) {
  vector<mpreal> data;

  if(numDigits < 1) numDigits = 1;
  if(numDigits > 15) numDigits = 15;

  mpreal precision = pow("10", -numDigits);

  string aa = "0.306892758965492";
  string bb = "0.125506782829762";

  string a1 = "0.";
  string b1 = "0.";
  for(int i=0; i<numDigits; i++) {
    a1 += aa[2+i];
    b1 += bb[2+i];
  }
  string a2 = "0.";
  string b2 = "0.";
  for(int i=0; i<numDigits+1; i++) {
    a2 += aa[2+i];
    b2 += bb[2+i];
  }

  mpreal a3 = atof(a1.c_str());
  mpreal b3 = atof(b1.c_str());

  mpreal a4 = atof(a2.c_str());
  mpreal b4 = atof(b2.c_str());

  mpreal diffa = a4-a3;
  mpreal diffb = b4-b3;

  if(diffa > precision/"2"-precision/"20") a3 += precision;
  if(diffb > precision/"2"-precision/"20") b3 += precision;

  data.push_back("1");
  data.push_back("-1");
  data.push_back("0");
  data.push_back("0");
  data.push_back(a3);
  data.push_back(b3);
  data.push_back("0");

  data.push_back("1");
  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back(a3);
  data.push_back(b3);
  data.push_back("0");

  data.push_back("1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("-2"*a3);
  data.push_back("-2"*b3);
  data.push_back("0");

  return data;
}

vector<mpreal> Initializer::get_pythagorean() {
  vector<mpreal> data;
  data.push_back("3");
  data.push_back("1");
  data.push_back("3");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back("4");
  data.push_back("-2");
  data.push_back("-1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back("5");
  data.push_back("1");
  data.push_back("-1");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  return data;
}
vector<mpreal> Initializer::get_sitnikov(mpreal ecc, mpreal vz0) {
  vector<mpreal> data = get_binary("1.0", ecc);
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back(vz0);
  return data;
}
vector<mpreal> Initializer::get_sitnikov(mpreal ecc, mpreal vz0, mpreal angle) {
  vector<mpreal> data = get_binary("1.0", ecc);
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back("0.0");
  data.push_back(vz0);

  rotate_x(data, angle);
  centralize(data);

  return data;
}

// N=4


// N=10
vector<mpreal> Initializer::get_solar_system() {
  vector<mpreal> data;
  data.push_back("1");
  data.push_back("-0.00712377699396443");
  data.push_back("-0.00283538985630487");
  data.push_back("-4.63675825294682e-06");
  data.push_back("0.000315043129753047");
  data.push_back("-0.000428788597909482");
  data.push_back("-8.53535576770256e-07");
      
  data.push_back("1.66013679527193e-07");
  data.push_back("-0.134041839403259");
  data.push_back("-0.450586209684483");
  data.push_back("-0.0317327781784651");
  data.push_back("1.24765113964416");
  data.push_back("-0.367348560236518");
  data.push_back("-0.11505457416379");
      
  data.push_back("2.44783833966455e-06");
  data.push_back("-0.726062021587661");
  data.push_back("-0.039665992657221");
  data.push_back("0.0218575042921458");
  data.push_back("0.0537323647367395");
  data.push_back("-1.17961128266637");
  data.push_back("-0.0273933407210185");
       
  data.push_back("3.04043264626853e-06");
  data.push_back("-0.189566509161968");
  data.push_back("0.963413224416489");
  data.push_back("0.00338899408041853");
  data.push_back("-0.998489725622231");
  data.push_back("-0.189356029048893");
  data.push_back("-0.0277905056358035");
      
  data.push_back("3.22715144505387e-07");
  data.push_back("1.38407082168322");
  data.push_back("-0.00854263708448956");
  data.push_back("0.00195338624778221");
  data.push_back("0.0339697602312997");
  data.push_back("0.882336954083445");
  data.push_back("0.0258975097250254");
      
  data.push_back("0.000954791938424327");
  data.push_back("3.97896476349658");
  data.push_back("2.95779610056192");
  data.push_back("0.0277878676318589");
  data.push_back("-0.267315308034173");
  data.push_back("0.372606717609675");
  data.push_back("0.000533607251937168");
      
  data.push_back("0.000285885980666131");
  data.push_back("6.37017339500169");
  data.push_back("6.60411921341949");
  data.push_back("-0.146003410253912");
  data.push_back("-0.2504985023359");
  data.push_back("0.224300955779407");
  data.push_back("0.00131013765082432");
      
  data.push_back("4.36625166911354e-05");
  data.push_back("14.502078723897");
  data.push_back("-13.6574724004636");
  data.push_back("0.0266710693728891");
  data.push_back("0.155097409556616");
  data.push_back("0.155754862016402");
  data.push_back("0.00394092485741823");
      
  data.push_back("5.15138902046612e-05");
  data.push_back("16.9343722652737");
  data.push_back("-24.904781203303");
  data.push_back("0.360719859921208");
  data.push_back("0.149752685968701");
  data.push_back("0.103693951061509");
  data.push_back("-0.000776459933629997");
      
  data.push_back("7.40740740830878e-09");
  data.push_back("-9.87695317832301");
  data.push_back("-28.0623724768707");
  data.push_back("5.35614463977912");
  data.push_back("0.177878490988979");
  data.push_back("-0.0885170840283217");
  data.push_back("-0.0375021230229398");
  
  // Normalize so that Msun=1, rpluto=1
  int N = data.size()/7;

  mpreal Msun = data[0];
  mpreal Cm = "1.0"/Msun;

  mpreal Rearth = sqrt(data[3*7+1]*data[3*7+1] + data[3*7+2]*data[3*7+2] + data[3*7+3]*data[3*7+3]);
  mpreal Cr = "1.0"/Rearth;

  mpreal Cv = sqrt(Cm/Cr);

  for(int i=0; i<N; i++) {
    data[i*7+0] *= Cm;
    data[i*7+1] *= Cr;
    data[i*7+2] *= Cr;
    data[i*7+3] *= Cr;
    data[i*7+4] *= Cv;
    data[i*7+5] *= Cv;
    data[i*7+6] *= Cv;
  }

  return data;
}

