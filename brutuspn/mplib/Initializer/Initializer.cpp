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
  if(config == "file_list_reversed") {  
    vector<mpreal> data;
    bool data_obtained = get_from_file_list_reversed(data, par[0], atoi(par[1].c_str()));
    if(data_obtained) return data;
  }
  if(config == "continue_from_file_reversed") {  
    vector<mpreal> data;
    bool data_obtained = continue_from_file_reversed(data, par[0]);
    if(data_obtained) return data;
  }
  if(config == "file_perturbed") {
    vector<mpreal> data;
    bool data_obtained = get_from_file_perturbed(data, par[0], atoi(par[1].c_str()), atoi(par[2].c_str()), par[3]);
    if(data_obtained) return data;
  }

  if(N == 2) {
    if(config == "circle") return get_circle( par[0] );
    if(config == "binary_random") return get_binary();
    if(config == "binary") {
      int M = par.size();
      if(M == 2) return get_binary(par[0], par[1]);
      if(M == 4) return get_binary(par[0], par[1], par[2], par[3] );
    }
    if(config == "parabola") return get_parabola(par[0]);
    if(config == "hyperbola") return get_hyperbola(par[0], par[1]);
    if(config == "radial") return get_radial(par[0], par[1]);
  }
  if(N == 3) {
    if(config == "figure8") return get_figure8();
    if(config == "figure8_precision") return get_figure8( atoi(par[0].c_str()) );
    if(config == "figure8_perturbed") return get_figure8_perturbed( atoi(par[0].c_str()), atoi(par[1].c_str()), par[2].c_str() );

    if(config == "triangle") return get_triangle();
    if(config == "triangle_precision") return get_triangle( atoi(par[0].c_str()) );
    if(config == "triangle_perturbed") return get_triangle_perturbed( atoi(par[0].c_str()), atoi(par[1].c_str()), par[2] );
    if(config == "triangle_scaled_perturbed") return get_triangle_scaled_perturbed( par[0], par[1], par[2], atoi(par[3].c_str()), atoi(par[4].c_str()), par[5] );

    if(config == "butterfly") return get_butterfly();
    if(config == "butterfly_precision") return get_butterfly( atoi(par[0].c_str()) );

    if(config == "pythagorean") return get_pythagorean();
    if(config == "pythagorean_perturbed") return get_pythagorean_perturbed( atoi(par[0].c_str()), atoi(par[1].c_str()), par[2].c_str() );

    if(config == "sitnikov") return get_sitnikov(par[0], par[1]);
    if(config == "sitnikov_angle") return get_sitnikov(par[0], par[1], par[2]);

    if(config == "sun_planet_comet") return get_sun_planet_comet_system(par[0]);
    if(config == "sun_planet_comet_perturbed") return get_sun_planet_comet_system_perturbed(par[0], atoi(par[1].c_str()));
    if(config == "sun_two_planets") return get_sun_two_planets_system(par[0], par[1]);
    if(config == "sun_two_planets_perturbed") return get_sun_two_planets_system_perturbed(par[0], par[1], atoi(par[2].c_str()));
    if(config == "agekyan_triple") return  get_agekyan_triple(par[0], par[1], atoi(par[2].c_str()), atoi(par[3].c_str()), atoi(par[4].c_str()) );
  }

  if(config == "exo_64443") {
    vector<mpreal> data;
    bool data_obtained = get_exo_64443(data, par[0]);
    if(data_obtained) return data;
  }
  if(config == "exo_64443_perturbed") {
    vector<mpreal> data;
    bool data_obtained = get_exo_64443_perturbed(data, par[0], atoi(par[1].c_str()), atoi(par[2].c_str()), par[3]);
    if(data_obtained) return data;
  }

  if(config == "solar_system_jpl_variable_N") {
    return get_solar_system_jpl_variable_N(atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()), atoi(par[3].c_str()), atoi(par[4].c_str()), atoi(par[5].c_str()), atoi(par[6].c_str()), atoi(par[7].c_str()), atoi(par[8].c_str()), atoi(par[9].c_str()), atoi(par[10].c_str()), atof(par[11].c_str()));
  }

  if(config == "plummer") {
    cerr << par.size() << endl;
    return get_plummer(  N, atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str())  );
  }
  if(config == "plummer_perturbed") {
    return get_plummer_perturbed(N, atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()), atoi(par[3].c_str()), atoi(par[4].c_str()), par[5]);
  }
  if(config == "plummer_static") {
    return get_static_plummer(N, atoi(par[0].c_str()), atoi(par[1].c_str()), atoi(par[2].c_str()), par[3], par[4], par[5], par[6]);
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
  int N = data.size()/7;
  if(N>1) {
    mpreal M = diagnostics.get_mass(data);
    mpreal Cm = "1"/M;
    for(int i=0; i<N; i++) {
      data[i*7] *= Cm;
    }

    mpreal r_vir = diagnostics.get_virial_radius(data);
    mpreal Cr = "1"/r_vir;
    for(int i=0; i<N; i++) {
      data[i*7+1] *= Cr;
      data[i*7+2] *= Cr;
      data[i*7+3] *= Cr;
    }

    mpreal sigma = diagnostics.get_velocity_dispersion(data);
    if(sigma > "0") {
      mpreal Cv = "1"/sqrt("2")/sigma;
      for(int i=0; i<N; i++) {
        data[i*7+4] *= Cv;
        data[i*7+5] *= Cv;
        data[i*7+6] *= Cv;
      }
    }   
  }
}
void Initializer::normalize_to_fraction_of_virial(vector<mpreal> &data, mpreal Q) {
  int N = data.size()/7;
  if(N>1) {
    mpreal ek = diagnostics.get_kinetic_energy(data);
    mpreal ep = diagnostics.get_potential_energy(data);

    mpreal Cv = sqrt(-Q*ep/ek);
    for(int i=0; i<N; i++) {
      data[i*7+4] *= Cv;
      data[i*7+5] *= Cv;
      data[i*7+6] *= Cv;
    }

    ek = diagnostics.get_kinetic_energy(data);
    ep = diagnostics.get_potential_energy(data);

    mpreal C = "-0.25" / (ek+ep);
    mpreal Cr = "1" / C;
    Cv = sqrt(C);

    for(int i=0; i<N; i++) {
      data[i*7+1] *= Cr;
      data[i*7+2] *= Cr;
      data[i*7+3] *= Cr;
      data[i*7+4] *= Cv;
      data[i*7+5] *= Cv;
      data[i*7+6] *= Cv;
    }   

    ek = diagnostics.get_kinetic_energy(data);
    ep = diagnostics.get_potential_energy(data);
    cerr << "e = " << ek+ep << endl;
    cerr << "q = " << ek/ep << endl;
  } 
}
void Initializer::normalize_to_boekholt_units(vector<mpreal> &data) {
  int N = data.size()/7;
  mpreal M = diagnostics.get_mass(data);
  mpreal Cm = "1"/M;
  for(int i=0; i<N; i++) {
    data[i*7] *= Cm;
  }
  mpreal r_har = diagnostics.get_harmonic_radius(data);
  mpreal Cr = "1"/r_har;
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
void Initializer::rescale_velocity(vector<mpreal> &data, mpreal Cv) {
  int N = data.size()/7;
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
  mpreal m_ratio = (mpreal)random.uniform(1.0, 10.0);
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
  cerr << data_handler.get_data().size() << " " << file << endl;
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[0];
    return true;
  }
}
bool Initializer::get_from_file_list(vector<mpreal> &data, string file, int index) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[index];
    return true;
  }
}
bool Initializer::continue_from_file(vector<mpreal> &data, string file) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[ file_data.size()-1 ];
    return true;
  }
}

bool Initializer::get_from_file_list_reversed(vector<mpreal> &data, string file, int index) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[index];
    flip_v(data);
    return true;
  }
}

bool Initializer::continue_from_file_reversed(vector<mpreal> &data, string file) {
  data.clear();
  bool file_read = data_handler.process(file);
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[ file_data.size()-1 ];
    flip_v(data);
    return true;
  }
}
bool Initializer::get_from_file_perturbed(vector<mpreal> &data, string file, int index, int coor, mpreal perturbation) {
  data.clear();
  bool file_read = data_handler.process(file);
  cerr << data_handler.get_data().size() << " " << file << endl;
  if(!file_read) return false;
  else {
    vector< vector<mpreal> > file_data = data_handler.get_data();
    data = file_data[0];
  
    data[index*7+coor] += perturbation;

    return true;
  }
}

// N=2
vector<mpreal> Initializer::get_circle(mpreal m_ratio) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1"/m_ratio;
  mpreal m2 = "1"/("1"+m_ratio);
  mpreal m1 = "1"-m2;

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

  if(m_ratio < "1") m_ratio = "1"/m_ratio;
  mpreal m2 = "1"/("1"+m_ratio);
  mpreal m1 = "1"-m2;

  mpreal ra = "1"+eccentricity;
  mpreal rp = "1"-eccentricity;

  mpreal v_aph = "1";
  if(eccentricity != "0") {
    v_aph = sqrt( "2" * ("1"/rp-"1"/ra) / (ra*ra/(rp*rp)-"1") );
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

  if(m_ratio < "1") m_ratio = "1"/m_ratio;
  mpreal m2 = "1"/("1"+m_ratio);
  mpreal m1 = "1"-m2;

  mpreal ra = "1"+eccentricity;
  mpreal rp = "1"-eccentricity;

  mpreal v_aph = "1";
  if(eccentricity != "0") {
    v_aph = sqrt( "2" * ("1"/rp-"1"/ra) / (ra*ra/(rp*rp)-"1") );
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
vector<mpreal> Initializer::get_binary(mpreal mu, mpreal m_ratio, mpreal semimajor_axis, mpreal eccentricity) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1"/m_ratio;
  mpreal m2 = "1"/("1"+m_ratio);
  mpreal m1 = "1"-m2;

  mpreal ra = "1"+eccentricity;
  mpreal rp = "1"-eccentricity;

  mpreal v_aph = "1";
  if(eccentricity != "0") {
    v_aph = sqrt( "2" * ("1"/rp-"1"/ra) / (ra*ra/(rp*rp)-"1") );
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

  mpreal Cm = mu;
  mpreal Cr = semimajor_axis;
  rescale_mass_position(data, Cm, Cr);  

  return data;
}

vector<mpreal> Initializer::get_parabola(mpreal m_ratio) {
  vector<mpreal> data;

  if(m_ratio < "1") m_ratio = "1"/m_ratio;
  mpreal m2 = "1"/("1"+m_ratio);
  mpreal m1 = "1"-m2;

  mpreal rp = "2";
  mpreal vp = "1";

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

  if(m_ratio > "1") m_ratio = "1"/m_ratio;
  mpreal m1 = "1"/("1"+m_ratio);
  mpreal m2 = "1"-m1;

  if(eccentricity <= "1") eccentricity = (mpreal)"1" + (mpreal)"1e-12";
  mpreal e = "0.5";

  mpreal vp = sqrt( "2"*(e + "1"/("1"-eccentricity)) );
  mpreal rp = "2" / (vp*vp - "1");

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

  if(m_ratio > "1") m_ratio = "1"/m_ratio;
  mpreal m1 = "1"/("1"+m_ratio);
  mpreal m2 = "1"-m1;

  mpreal x1 = "0";
  mpreal x2 = "1e-6";

  mpreal vx1 = "0";
  mpreal vx2 = sqrt("2"*(energy+"1"/x2));

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

  mpreal precision = pow("10.0", (mpreal)-numDigits);

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
vector<mpreal> Initializer::get_figure8_perturbed(int index, int coor, mpreal delta) {
  vector<mpreal> data = get_figure8();
  data[index*7 + coor] += delta;
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
  data.push_back(-sqrt("0.75"));
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

  mpreal precision = pow("10.0", (mpreal)-numDigits);

  string cc = "0.8660254037844386";

  string c1 = "0.";
  for(int i=0; i<numDigits; i++) {
    c1 += cc[2+i];
  }
  string c2 = "0.";
  for(int i=0; i<numDigits+1; i++) {
    c2 += cc[2+i];
  }

  mpreal c3 = c1;
  mpreal c4 = c2;

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
vector<mpreal> Initializer::get_triangle_perturbed(int index, int coor, mpreal dx) {
  vector<mpreal> data;
  data = get_triangle();  
  data[index*7+coor] += dx;
  return data;
}

vector<mpreal> Initializer::get_triangle_scaled_perturbed(mpreal Cm, mpreal Cr, mpreal Cv, int index, int coor, mpreal dx) {
  vector<mpreal> data;
  data = get_triangle();
  
  //mpreal Cv = sqrt(Cm/Cr);
  for(int i=0; i<3; i++) {
    data[i*7+0] *= Cm;
    data[i*7+1] *= Cr;
    data[i*7+2] *= Cr;
    data[i*7+3] *= Cr;
    data[i*7+4] *= Cv;
    data[i*7+5] *= Cv;
    data[i*7+6] *= Cv;
  }

  data[index*7+coor] += dx;

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

  mpreal precision = pow("10.0", (mpreal)-numDigits);

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

  mpreal a3 = a1;
  mpreal b3 = b1;

  mpreal a4 = a2;
  mpreal b4 = b2;

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
vector<mpreal> Initializer::get_pythagorean_perturbed(int index, int coor, mpreal delta) {
  vector<mpreal> data = get_pythagorean();
  data[index*7 + coor] += delta;
  return data;
}

vector<mpreal> Initializer::get_sitnikov(mpreal ecc, mpreal vz0) {
  vector<mpreal> data = get_binary((mpreal)"1", ecc);
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
  vector<mpreal> data = get_binary((mpreal)"1", ecc);
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

vector<mpreal> Initializer::get_agekyan_triple(mpreal x, mpreal y, int seed, int pivot, int index) {
  set_random_generator(seed, pivot);

  vector<mpreal> data;

  data.push_back("1");
  data.push_back("-0.5");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  data.push_back("1");
  data.push_back("0.5");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");
  data.push_back("0");

  if(x < "0" || y < "0") {
    mpreal x0 = "0";
    mpreal y0 = "0";
    mpreal z0 = "0";

    for(int i=0; i<index+1; i++) {
      while (1) {
        x0 = (mpreal)random.uniform(0., 0.5);
        y0 = (mpreal)random.uniform(0., 0.86602540378);
        z0 = (x0 + 0.5)*(x0 + 0.5) + y0*y0;
        if(z0 < "1") {
          break;        
        }
      }
    }

    data.push_back("1");
    data.push_back(x0);
    data.push_back(y0);
    data.push_back("0");
    data.push_back("0");
    data.push_back("0");
    data.push_back("0");
  }
  else {
    data.push_back("1");
    data.push_back(x);
    data.push_back(y);
    data.push_back("0");
    data.push_back("0");
    data.push_back("0");
    data.push_back("0");
  }
  cerr << data[2*7+1] << " " << data[2*7+2] << endl;
  return data;
}

vector<mpreal> Initializer::get_sun_planet_comet_system(mpreal MA) {
  mpreal a2_1 = (mpreal)pow((mpreal)"2", (mpreal)"2"/"3");
  mpreal pi = acos("-1");

  mpreal M   = "1";
  mpreal x0  = "0";
  mpreal y0  = "0";
  mpreal z0  = "0";
  mpreal vx0 = "0";
  mpreal vy0 = "0";
  mpreal vz0 = "0";

  mpreal m1  = "1.e-3";
  mpreal x1  = "1";
  mpreal y1  = "0";
  mpreal z1  = "0";
  mpreal vx1 = "0";
  mpreal vy1 = sqrt((M+m1)/x1);
  mpreal vz1 = "0";

  mpreal m2  = "1.e-16";
  mpreal x2  = a2_1*cos(MA/180*pi);
  mpreal y2  = a2_1*sin(MA/180*pi);
  mpreal z2  = "0";
  mpreal vx2 = sqrt((M+m2)/a2_1)*-sin(MA/180*pi);
  mpreal vy2 = sqrt((M+m2)/a2_1)*cos(MA/180*pi);
  mpreal vz2 = "0";

  vector<mpreal> data;

  data.push_back(M);
  data.push_back(x0);
  data.push_back(y0);
  data.push_back(z0);
  data.push_back(vx0);
  data.push_back(vy0);
  data.push_back(vz0);

  data.push_back(m1);
  data.push_back(x1);
  data.push_back(y1);
  data.push_back(z1);
  data.push_back(vx1);
  data.push_back(vy1);
  data.push_back(vz1);

  data.push_back(m2);
  data.push_back(x2);
  data.push_back(y2);
  data.push_back(z2);
  data.push_back(vx2);
  data.push_back(vy2);
  data.push_back(vz2);

  centralize(data);

  return data;
}

vector<mpreal> Initializer::get_sun_planet_comet_system_perturbed(mpreal MA, int index) {
  vector<mpreal> data = get_sun_planet_comet_system(MA);

  set_random_generator(0, 0);
  mpreal perturbation_size = "1e-14";

  mpreal dx = "0";
  mpreal dy = "0";
  mpreal dz = "0";
  mpreal dvx = "0";
  mpreal dvy = "0";
  mpreal dvz = "0";
  for(int i=0; i<index+1; i++) {
    dx  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dy  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dz  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvx = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvy = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvz = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
  }

  data[2*7+1] += dx;
  data[2*7+2] += dy;
  data[2*7+3] += dz;
  data[2*7+4] += dvx;
  data[2*7+5] += dvy;
  data[2*7+6] += dvz;

  return data;
}

vector<mpreal> Initializer::get_sun_two_planets_system(mpreal m, mpreal MA) {
  mpreal a2_1 = (mpreal)pow((mpreal)"2", (mpreal)"2"/"3");
  mpreal pi = acos("-1");

  mpreal M   = "1";
  mpreal x0  = "0";
  mpreal y0  = "0";
  mpreal z0  = "0";
  mpreal vx0 = "0";
  mpreal vy0 = "0";
  mpreal vz0 = "0";

  mpreal m1  = m;
  mpreal x1  = "1";
  mpreal y1  = "0";
  mpreal z1  = "0";
  mpreal vx1 = "0";
  mpreal vy1 = sqrt((M+m1)/x1);
  mpreal vz1 = "0";

  mpreal m2  = m;
  mpreal x2  = a2_1*cos(MA/180*pi);
  mpreal y2  = a2_1*sin(MA/180*pi);
  mpreal z2  = "0";
  mpreal vx2 = sqrt((M+m2)/a2_1)*-sin(MA/180*pi);
  mpreal vy2 = sqrt((M+m2)/a2_1)*cos(MA/180*pi);
  mpreal vz2 = "0";

  vector<mpreal> data;

  data.push_back(M);
  data.push_back(x0);
  data.push_back(y0);
  data.push_back(z0);
  data.push_back(vx0);
  data.push_back(vy0);
  data.push_back(vz0);

  data.push_back(m1);
  data.push_back(x1);
  data.push_back(y1);
  data.push_back(z1);
  data.push_back(vx1);
  data.push_back(vy1);
  data.push_back(vz1);

  data.push_back(m2);
  data.push_back(x2);
  data.push_back(y2);
  data.push_back(z2);
  data.push_back(vx2);
  data.push_back(vy2);
  data.push_back(vz2);

  centralize(data);

  return data;

}
vector<mpreal> Initializer::get_sun_two_planets_system_perturbed(mpreal m, mpreal MA, int index) {
  vector<mpreal> data = get_sun_two_planets_system(m, MA);

  set_random_generator(0, 0);
  mpreal perturbation_size = "1e-14";

  mpreal dx = "0";
  mpreal dy = "0";
  mpreal dz = "0";
  mpreal dvx = "0";
  mpreal dvy = "0";
  mpreal dvz = "0";
  for(int i=0; i<index+1; i++) {
    dx  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dy  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dz  = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvx = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvy = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
    dvz = (mpreal)random.uniform(-1.0, 1.0) * perturbation_size;
  }

  data[2*7+1] += dx;
  data[2*7+2] += dy;
  data[2*7+3] += dz;
  data[2*7+4] += dvx;
  data[2*7+5] += dvy;
  data[2*7+6] += dvz;

  return data;
}

vector<mpreal> Initializer::get_solar_system_jpl_variable_N(bool Mercury, bool Venus, bool Earth, bool Mars, bool Jupiter, bool Saturn, bool Uranus, bool Neptune, bool Pluto, int index, int coor, mpreal perturbation) {
  vector<mpreal> data;

  data.push_back("1");
  data.push_back("2.072946918590417E-03");
  data.push_back("-1.643061287382110E-03");
  data.push_back("-1.185823088540438E-04");
  data.push_back("5.059930816518583E-06");
  data.push_back("3.989470642047999E-06");
  data.push_back("-1.243838264183300E-07");

  if(Mercury) {      
    data.push_back("1.6601375118291464e-7");
    data.push_back("3.574733791200970E-01");
    data.push_back("-1.171187317227507E-02");
    data.push_back("-3.354856363159710E-02");
    data.push_back("-4.633283577618537E-03");
    data.push_back("2.938619087935777E-02");
    data.push_back("2.826185769428105E-03");
  }
  if(Venus) {      	   	       
    data.push_back("2.44783857411814e-6");
    data.push_back("4.876790422110325E-01");
    data.push_back("5.332660900104316E-01");
    data.push_back("-2.081177373810901E-02");
    data.push_back("-1.503296302925791E-02");
    data.push_back("1.351327482426259E-02");
    data.push_back("1.052862263953481E-03");
  }
  if(Earth) {       	  	         
    data.push_back("3.0404326668433355e-6");
    data.push_back("4.026754002834322E-01");
    data.push_back("-9.357700005133348E-01");
    data.push_back("-8.932490601085987E-05");
    data.push_back("1.553740564985712E-02");
    data.push_back("6.720378960900881E-03");
    data.push_back("-3.796916120995478E-07");
  }	  	       
  if(Mars) {      
    data.push_back("3.227150370113422e-7");
    data.push_back("-5.844200241438797E-01");
    data.push_back("-1.390641490399950E+00");
    data.push_back("-1.482626102558087E-02");
    data.push_back("1.342616848375762E-02");
    data.push_back("-4.237723031036153E-03");
    data.push_back("-4.184228274540421E-04");
  }
  if(Jupiter) {      
    data.push_back("0.0009545940905009727");
    data.push_back("-2.705363190538412E+00");
    data.push_back("4.511969711234348E+00");
    data.push_back("4.171926801427295E-02");
    data.push_back("-6.562761960592283E-03");
    data.push_back("-3.524082083761892E-03");
    data.push_back("1.615059047295873E-04");
  }
  if(Saturn) {      	   	       
    data.push_back("0.0002858150131819827");
    data.push_back("-6.120971157013838E+00");
    data.push_back("-7.799968424277239E+00");
    data.push_back("3.792260890576647E-01");
    data.push_back("4.084918553333202E-03");
    data.push_back("-3.458702472764618E-03");
    data.push_back("-1.026661638527525E-04");
  }
  if(Uranus) {        	   	        
    data.push_back("4.3658047415668107e-5");
    data.push_back("1.947563780483501E+01");
    data.push_back("4.637428187478366E+00");
    data.push_back("-2.350901336659057E-01");
    data.push_back("-9.397988346434892E-04");
    data.push_back("3.642810279796556E-03");
    data.push_back("2.574368503661651E-05");
  }
  if(Neptune) {      	  	       
    data.push_back("5.1503137142531555e-5");
    data.push_back("2.731924333900583E+01");
    data.push_back("-1.233234986862218E+01");
    data.push_back("-3.756379750691105E-01");
    data.push_back("1.270955946848631E-03");
    data.push_back("2.880167515810431E-03");
    data.push_back("-8.840273274466099E-05");
  }
  if(Pluto) {          	   	   
    data.push_back("7.3442323846354385e-9");
    data.push_back("6.872419409625540E+00");
    data.push_back("-3.192773277605658E+01");
    data.push_back("1.428555846904559E+00");
    data.push_back("3.130439755573468E-03");
    data.push_back("2.364247270061568E-05");
    data.push_back("-9.080391301223528E-04");
  }

  //mpreal Cv = 58.12435834603945838;
  mpreal Cv = "365.242" / ("2"*acos("-1."));

  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=4; j<7; j++) {
      data[i*7+j] *= Cv;
    }
  }

  centralize(data);

  data[index*7+coor] += perturbation;

  return data;
}

// Exoplanets
bool Initializer::get_exo_64443(vector<mpreal> &data, string file) {
    bool loaded = get_from_file(data, file);
    if(loaded == true) {
        mpreal Cv = "365.25" / ("2"*acos("-1")); //"2"*acos("-1"); // x 2*pi*AU/yr = C*y AU/d, C = 365.25  
        rescale_velocity(data, Cv);
    }
    return loaded;
}
bool Initializer::get_exo_64443_perturbed(vector<mpreal> &data, string file, int particle, int coor, mpreal delta) {
    bool loaded = get_exo_64443(data, file);
    if(loaded == true) {
        data[particle*7+coor] += delta;
    }
    return loaded;
}

vector<mpreal> Initializer::get_plummer(int N, int seed, int pivot, int index) {
  vector<mpreal> data(N*7);

  set_random_generator(seed, pivot);

  for(int i=0; i<=index; i++) {
    vector<mpreal> m(N), x(N), y(N), z(N), vx(N), vy(N), vz(N);
    mpreal pi = acos("-1");

    // Set masses
    mpreal c = N;
    for(int i=0; i<N; i++) {
      m[i] = "1"/c;
    } 

    // Set positions
    for(int i=0; i<N; i++) {
      mpreal a = "-2";
      mpreal b = "3";
      mpreal r = "1" / sqrt( pow((mpreal)random.uniform(0.0, 1.0), a/b ) - "1.0");
      mpreal theta = acos( (mpreal)random.uniform(0.0, 1.0)*"2.0" - "1.0");
      mpreal phi = (mpreal)random.uniform(0.0, 1.0) * "2.0"*pi;
      x[i] = r*sin(theta)*cos(phi);
      y[i] = r*sin(theta)*sin(phi);
      z[i] = r*cos(theta);
    } 

    // Set velocities
    for(int i=0; i<N; i++) {    
      mpreal xx = "0.0";
      mpreal yy = "0.1";
      while( yy > xx*xx*pow( ("1.0"-xx*xx), "3.5" ) ) {
        xx = (mpreal)random.uniform(0.0, 1.0);
        yy = (mpreal)random.uniform(0.0, 1.0) / "10.0";
      }
      mpreal r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
      mpreal v = xx * sqrt("2.0") * pow( ( "1.0" + r2), "-0.25" );

      mpreal theta = acos( (mpreal)random.uniform(0.0, 1.0) * "2.0" - "1.0" );
      mpreal phi = (mpreal)random.uniform(0.0, 1.0) * "2.0"*pi;

      vx[i] = v*sin(theta)*cos(phi);
      vy[i] = v*sin(theta)*sin(phi);
      vz[i] = v*cos(theta);
    }

    for(int i=0; i<N; i++) {
      data[i*7] = m[i];
      data[i*7+1] = x[i];
      data[i*7+2] = y[i];
      data[i*7+3] = z[i];
      data[i*7+4] = vx[i];
      data[i*7+5] = vy[i];
      data[i*7+6] = vz[i];
    }
  } 

  centralize(data);
  normalize_to_nbody_units(data);

  return data;
}
vector<mpreal> Initializer::get_plummer_perturbed(int N, int seed, int pivot, int index, int star, int coor, mpreal delta) {
  vector<mpreal> data = get_plummer(N, seed, pivot, index);

  data[star*7+coor] += delta;

  return data;
}

vector<mpreal> Initializer::get_static_plummer(int N, int seed, int pivot, int index, mpreal m0, mpreal x0, mpreal vz0, mpreal delta) {
  //vector<mpreal> orbit = get_plummer(1000, seed+1, pivot, index);
  vector<mpreal> bodies = get_plummer(N, seed, pivot, index);

  vector<mpreal> stars;

  mpreal m = m0;
  mpreal x = x0;
  mpreal y = "0";
  mpreal z = "0";
  mpreal vx = "0";
  mpreal vy = "0";
  mpreal vz = vz0;

  stars.push_back(m);
  stars.push_back(x);
  stars.push_back(y);
  stars.push_back(z);
  stars.push_back(vx);
  stars.push_back(vy);
  stars.push_back(vz);

  //for(int i=orbitIndex*7; i<(orbitIndex+1)*7; i++) {
  //  stars.push_back(orbit[i]);
  //}
  for(int i=0; i<bodies.size(); i++) {
    stars.push_back(bodies[i]);
  }

  //stars[0] = m0;

  int coor = 1;
  stars[coor] += delta;

  return stars;
}








