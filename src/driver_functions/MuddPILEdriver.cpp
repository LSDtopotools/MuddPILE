//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// MuddPILEdriver
// The main driver function for the MuddPILE model
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Copyright (C) 2017 Simon M. Mudd 2017
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sys/stat.h>
#include "../LSDRasterModel.hpp"
#include "../LSDParticleColumn.hpp"
#include "../LSDParticle.hpp"
#include "../LSDCRNParameters.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDStatsTools.hpp"
#include "../LSDRasterMaker.hpp"
using namespace std;


int main (int nNumberofArgs,char *argv[])
{

  //Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "================================================================" << endl;
    cout << "|| Welcome to the MuddPILEdirver!                             ||" << endl;
    cout << "|| This program drives the MuddPILE lanscape evolution model. ||" << endl;
    cout << "|| One day this model will have documentation.                ||" << endl;
    cout << "|| This program was developed by Simon M. Mudd                ||" << endl;
    cout << "||  at the University of Edinburgh                            ||" << endl;
    cout << "================================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be: " << endl;
    cout << "In linux:" << endl;
    cout << "./MuddPILEdriver.exe /LSDTopoTools/Topographic_projects/Model_runs/ MuddPILE_test.driver" << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }

  string path_name = argv[1];
  string f_name = argv[2];

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // Parameters for initiating the model domain
  int_default_map["NRows"] = 200;
  int_default_map["NCols"] = 400;
  float_default_map["DataResolution"] = 30;

  // Parameters that are used if you load a raster
  bool_default_map["read_initial_raster"] = false;
  float_default_map["minimum_elevation"] = 0.0;
  float_default_map["maximum_elevation"] = 30000;
  float_default_map["min_slope_for_fill"] = 0.0001;
  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  bool_default_map["print_raster_without_seas"] = false;

  // paramters for controlling model output
  int_default_map["print_interval"] = 10;
  bool_default_map["write_hillshade"] = true;

  // control of m and n, and paramters for chi
  float_default_map["A_0"] = 1;
  float_default_map["m"] = 0.5;
  float_default_map["n"] = 1;
  int_default_map["uplift_mode"] = 0;
  float_default_map["dt"] = 250;
  float_default_map["D"] = 0.002;
  float_default_map["S_c"] = 1.0;
  float_default_map["background_K"] = 0.0005;

  // Parameters for the initial surface
  bool_default_map["use_diamond_square_initial"] = true;
  float_default_map["diamond_square_relief"] = 16;
  int_default_map["diamond_square_feature_order"] = 8;
  int_default_map["ds_minimum_feature_order"] = 4;
  bool_default_map["taper_edges"] = true;
  int_default_map["taper_rows"] = 10;
  bool_default_map["superimpose_parabola"] = true;
  float_default_map["parabola_relief"] = 6;
  bool_default_map["roughen_surface"] = true;
  bool_default_map["fill_then_roughen_surface"] = true;
  float_default_map["roughness_relief"] = 0.25;
  bool_default_map["diffuse_initial_surface"] = false;
  int_default_map["diffuse_steps"] = 5;
  bool_default_map["print_initial_surface"] = true;
  

  // Parameters for spinning up the simulation
  bool_default_map["spinup"] = false;
  float_default_map["spinup_K"] = 0.001;
  float_default_map["spinup_U"] = 0.001;
  float_default_map["spinup_dt"] = 250;
  float_default_map["spinup_time"] = 20000;
  bool_default_map["staged_spinup"] = true;
  
  // This spinup routine tries to combine snapping and hillslope diffusion
  bool_default_map["cyclic_spinup"] = false;
  int_default_map["spinup_cycles"] = 5;
  bool_default_map["cycle_K"] = true;
  bool_default_map["cycle_U"] = false;
  float_default_map["cycle_K_factor"] = 2;
  float_default_map["cycle_U_factor"] = 2;
  
  // The diamond square spinup routine
  // This generates the nicest initial surfaces. 
  bool_default_map["diamond_square_spinup"] = false;

  // control of snapping to steady state
  bool_default_map["snap_to_steady"] = false;
  float_default_map["snapped_to_steady_uplift"] = 0.0001;
  float_default_map["snapped_to_steady_relief"] = 400;
  bool_default_map["print_snapped_to_steady_frame"] = false;
  
  // forcing of dissection
  bool_default_map["force_dissect"] = true;
  int_default_map["force_dissect_steps"] = 10000;

  // Some parameters for very rudimentary steady forcing
  bool_default_map["rudimentary_steady_forcing"] = false;
  float_default_map["rudimentary_steady_forcing_time"] = 100000;
  float_default_map["rudimentary_steady_forcing_uplift"] = 0.0005;
  float_default_map["rudimentary_steady_forcing_K"] = 0.0005;

  // Parameters for hillslopes
  bool_default_map["hillslopes_on"] = false;

  // Some parameters for cyclic forcing
  // these also inherit parameters from the cyclic spinup
  bool_default_map["run_cyclic_forcing"] = false;
  float_default_map["cyclic_forcing_time"] = 10000;
  float_default_map["baseline_K_for_cyclic"] = 0.0001;
  float_default_map["baseline_U_for_cyclic"] = 0.0005;
  int_default_map["cyclic_cycles"] = 3;
  float_default_map["cyclic_dt"] = 250;
  

  // some parameters for setting K to a fixed uplift and relief
  bool_default_map["set_fixed_relief"] = false;
  float_default_map["fixed_relief"] = 500;
  
  
  // some parameters for having a random uplift through time
  bool_default_map["run_random_forcing"] = false;
  float_default_map["maximum_time_for_random_cycle"] = 20000;
  float_default_map["minimum_time_for_random_cycle"] = 5000;
  float_default_map["maximum_U_for_random_cycle"] = 0.001;
  float_default_map["minimum_U_for_random_cycle"] = 0.0001;
  float_default_map["random_dt"] = 10;  
  int_default_map["random_cycles"] = 4;
  
  // some parameters for a spatially varying K and U
  bool_default_map["make_spatially_varying_K"] = false;
  float_default_map["spatially_varying_max_K"] = 0.0001;
  float_default_map["spatially_varying_min_K"] = 0.000001;
  int_default_map["min_blob_size"] = 50;
  int_default_map["max_blob_size"] = 100;
  int_default_map["n_blobs"] = 10;
  
  bool_default_map["spatially_varying_forcing"] = false;
  bool_default_map["spatially_varying_K"] = false;
  bool_default_map["spatially_varying_U"] = false;
  bool_default_map["calculate_K_from_relief"] = false;
  int_default_map["spatial_K_method"] = 0;
  int_default_map["spatial_U_method"] = 1;          // 0 doesn't work
  bool_default_map["load_K_raster"] = false;
  bool_default_map["load_U_raster"] = false;
  
  float_default_map["spatial_K_factor"] = 3;
  float_default_map["spatial_variation_time"] = 20000;
  float_default_map["min_U_for_spatial_var"] = 0.0001;
  float_default_map["max_U_for_spatial_var"] = 0.0005;
  int_default_map["K_smoothing_steps"] = 2;
  float_default_map["spatial_dt"] = 100;
  int_default_map["spatial_cycles"] = 5;
  bool_default_map["use_adaptive_timestep"] = true;
  float_default_map["maximum_timestep"] = 500;
  float_default_map["float_print_interval"] = 2000;
  bool_default_map["snap_to_steep_for_spatial_uplift"] = false;
  bool_default_map["snap_to_minimum_uplift"] = true;
  

  // Use the parameter parser to get the maps of the parameters required for the analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  // Now print the parameters for bug checking
  LSDPP.print_parameters();

  // location of the files
  string DATA_DIR =  LSDPP.get_read_path();
  string DEM_ID =  LSDPP.get_read_fname();
  string OUT_DIR = LSDPP.get_write_path();
  string OUT_ID = LSDPP.get_write_fname();
  string raster_ext =  LSDPP.get_dem_read_extension();
  vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
  string CHeads_file = LSDPP.get_CHeads_file();

  // Initiate the model object
  LSDRasterModel mod;

  cout << "Write filename is: " << OUT_DIR+OUT_ID << endl;
  mod.add_path_to_names(OUT_DIR);
  mod.set_name(OUT_DIR+OUT_ID);

  // set the print intervals and other parameters
  mod.set_m(this_float_map["m"]);
  mod.set_n(this_float_map["n"]);
  mod.set_uplift_mode(this_int_map["uplift_mode"]);
  mod.set_print_hillshade(this_bool_map["write_hillshade"]);
  mod.set_timeStep( this_float_map["dt"] );
  mod.set_print_interval(this_int_map["print_interval"]);
  mod.set_D( this_float_map["D"]);
  mod.set_S_c( this_float_map["S_c"] );
  mod.set_K( this_float_map["background_K"]);
  
  // print parameters to screen
  mod.print_parameters();

  // need this to keep track of the end time
  float current_end_time = 0;

  //============================================================================
  // Logic for reading an initial surface
  // It will look for a raster but if it doesn't find one it will default back to
  // fixed rows and columns. If you do load a raster it will override any
  // instructions about NRows, NCols, and DataResolution
  //============================================================================
  bool create_initial_surface = false;
  if(this_bool_map["read_initial_raster"])
  {
    cout << "I am going to try to read an intial raster for you." << endl;
    cout << "The read filename is: " <<  DATA_DIR+DEM_ID << endl;
    cout << "I am going to IGNORE initial surface instructions!" << endl;
    string header =  DATA_DIR+DEM_ID+".hdr";
    cout << "The full read path is: " << header << endl;
    ifstream file_info_in;
    file_info_in.open(header.c_str());
    // check if the parameter file exists
    if( not file_info_in.fail() )
    {
      cout << "I found the header. I am loading this initial file. " << endl;
      LSDRaster temp_raster(DATA_DIR+DEM_ID,"bil");
      if (this_bool_map["remove_seas"])
      {
        cout << "I am removing high and low values to get rid of things that should be nodata." << endl;
        float lower_threshold = this_float_map["minimum_elevation"];
        float upper_threshold = this_float_map["maximum_elevation"];
        bool belowthresholdisnodata = true;
        LSDRaster Flooded = temp_raster.mask_to_nodata_using_threshold(lower_threshold,belowthresholdisnodata);
        belowthresholdisnodata = false;
        temp_raster = Flooded.mask_to_nodata_using_threshold(upper_threshold,belowthresholdisnodata);
        if (this_bool_map["print_raster_without_seas"])
        {
          cout << "I'm replacing your input raster with a raster without seas." << endl;
          string this_raster_name = OUT_DIR+OUT_ID;
          temp_raster.write_raster(this_raster_name,raster_ext);
        }
      }

      LSDRasterModel temp_mod(temp_raster);
      mod = temp_mod;
    }
    else
    {
      cout << "Header doesn't exist. I am switching to a fixed rows and columns." << endl;
      cout << "I will also create an initial surface for you." << endl;
      create_initial_surface = true;
    }
  }
  else
  {
    cout << "You have chosen not to read an initial raster so I will create an initial surface for you." << endl;
    create_initial_surface = true;
  }


  // a bit of logic to override spinup logic if you choose the diamond square spinup
  if(this_bool_map["diamond_square_spinup"])
  {
    cout << "You have chosen the diamond square spinup. This orverrides all other spinup otions" << endl;
    cout << "It also overrides initial surface options" << endl;
    create_initial_surface = false;
    this_bool_map["cyclic_spinup"] = false;
    this_bool_map["spinup"] = false;
  }


  //============================================================================
  // Logic for creating an initial surface
  // There are a number of options here, but the defaults are the ones we have
  // found to create the nicest initial surfaces.
  // The diamond square routine creates a pseudo-fractal surface. The algorithm
  // comes from Minecraft's creator, believe it or not.
  // If you use only the diamond square, you will get pits in the middle of the
  //  intial surface which will fill, leaving areas of the DEM with very
  //  straight channels. To mitigate this you can superimpose a parabola
  //  on the surface, and you can also run a fill function and then add roughness
  //  to the filled DEM (these are the defaults).
  // The greater the relief of the parabola relative to the relief of the
  //  diamond square fractal, the straighter your channels and basins will be
  //  after dissection.
  //============================================================================
  if(create_initial_surface)
  {
    cout << endl << endl << "=============================================" << endl;
    cout << "Creating and initial surface. " << endl;
    cout << "NRows: " << this_int_map["NRows"]<< ", NCols: " << this_int_map["NCols"]
         << ", DataResolution: " << this_float_map["DataResolution"] << endl;
    mod.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"]);

    if (this_bool_map["use_diamond_square_initial"])
    {
      cout << "   I am starting with a fractal surface created by the diamond square algorithm." << endl;
      cout << "   It has a relief of " << this_int_map["diamond_square_relief"] << endl;
      cout << "   Let me check on the largest possible scale of features in the pseudo fractal surface." << endl;

      // we need to check the feature order to make sure it is not bigger than the largest dimension of the raster
      int smallest_dimension;
      if( this_int_map["NRows"] <= this_int_map["NCols"])
      {
        smallest_dimension = this_int_map["NRows"];
      }
      else
      {
        smallest_dimension = this_int_map["NCols"];
      }
      // now get the largest possible dimension
      float lps = floor( log2( float(smallest_dimension) ) );
      int largest_possible_scale = int(lps);
      if( this_int_map["diamond_square_feature_order"] > largest_possible_scale )
      {
        this_int_map["diamond_square_feature_order"] = largest_possible_scale;
      }
      cout << "    The largest dimnesion is: " << smallest_dimension << " pixels." << endl;
      cout << "    So the biggest scale is: " << pow(2,largest_possible_scale) << " pixels or 2 to the " << largest_possible_scale << " power" << endl;
      cout << "    I am setting the scale to: " << this_int_map["diamond_square_feature_order"] << endl;

      // Now actually build the fractal surface
      mod.intialise_diamond_square_fractal_surface(this_int_map["diamond_square_feature_order"], this_float_map["diamond_square_relief"]);

      if (this_bool_map["taper_edges"])
      {
        cout << "   Let me taper the edges of this fractal surface for you so that everything drains to the edge." << endl;
        cout << "   I am tapering along the " << this_int_map["taper_rows"] << " row closes to the N and S boundaries" << endl;
        mod.initialise_taper_edges_and_raise_raster(this_int_map["taper_rows"]);
      }
    }
    if (this_bool_map["superimpose_parabola"])
    {
      cout << "   I am superimposing a parabola with a relief of "  << this_float_map["parabola_relief"] << " metres" << endl;
      mod.superimpose_parabolic_surface(this_float_map["parabola_relief"]);
    }
    if(this_bool_map["roughen_surface"])
    {
      cout << "   I am going to roughen the surface for you. Roughness elements will " << endl;
      cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
      mod.set_noise(this_float_map["roughness_relief"]);
      mod.random_surface_noise();
    }
    if(this_bool_map["fill_then_roughen_surface"])
    {
      cout << "   I am going to fill and then roughen the surface for you. Roughness elements will " << endl;
      cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
      mod.raise_and_fill_raster();
      mod.set_noise(this_float_map["roughness_relief"]);
      mod.random_surface_noise();
    }
    
    if(this_bool_map["diffuse_initial_surface"])
    {
      cout << "I am going to diffuse the initial surface for you." << endl;
      mod.set_timeStep( 0.5 );
      for (int i = 0; i< this_int_map["diffuse_steps"]; i++)
      {
        cout << "Diffuse step " << i <<"/" << this_int_map["diffuse_steps"] << endl;
        mod.MuddPILE_nl_soil_diffusion_nouplift();
      }
      mod.set_timeStep( this_float_map["dt"] );
      mod.raise_and_fill_raster();
    }
    
    if(this_bool_map["print_initial_surface"])
    {
      int this_frame = 9999;
      mod.print_rasters_and_csv( this_frame );
    
    }
    cout << "Finished creating an initial surface. " << endl;
    cout << "=============================================" << endl << endl << endl;
  }

  //============================================================================
  // Logic for making spatially varying K fields
  //============================================================================
  if(this_bool_map["make_spatially_varying_K"])
  {
    cout << "I am going to make a spatially varying raster for the K parameter." << endl;
    
    LSDRasterMaker KRaster(this_int_map["NRows"],this_int_map["NCols"]);
    KRaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["spatially_varying_min_K"]);
    KRaster.random_square_blobs(this_int_map["min_blob_size"], this_int_map["max_blob_size"], 
                                this_float_map["spatially_varying_min_K"], this_float_map["spatially_varying_max_K"],
                                this_int_map["n_blobs"]);
  
  // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_KRaster";
    string bil_name = "bil";
    
    KRaster.write_raster(K_fname,bil_name);
    
    
    // now do a sine version. 
    LSDRasterMaker URaster(this_int_map["NRows"],this_int_map["NCols"]);
    URaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["spatially_varying_min_K"]);
    vector<float> x_coeff;
    vector<float> y_coeff;
    y_coeff.push_back(10);
    //y_coeff.push_back(0);
    //y_coeff.push_back(2);
    //y_coeff.push_back(0);
    //y_coeff.push_back(5);
    
    URaster.sine_waves(x_coeff, y_coeff);
    
    string U_fname = OUT_DIR+OUT_ID+"_URaster";
    URaster.write_raster(U_fname,bil_name);
  }



  //============================================================================
  // Logic for a spinning up the model.
  // Here spinup is used to dissect the landscape. Usually you create a
  // fractal surface first and then dissect it to get an initial channel and
  // basin geometry. You can then snap it to steady state with the snapping
  // functions.
  // This uses fluvial only! You need to run hillslope diffusion after this
  // if you want a steady landscape with hillslopes
  //============================================================================
  if(this_bool_map["spinup"])
  {
    cout << "I am going to spin the model up for you." << endl;
    cout << "This will rapidly develop a drainage network. It uses only fluvial incision." << endl;
    cout << "The typical spinup method is to dissect the landscape and then bring it to " << endl;
    cout << "steady state using the snap_to_steady flag." << endl;

    // turn the hillslope diffusion off
    mod.set_hillslope(false);

    current_end_time = this_float_map["spinup_time"];

    mod.set_endTime(this_float_map["spinup_time"]);
    mod.set_timeStep( this_float_map["spinup_dt"] );
    mod.set_K(this_float_map["spinup_K"]);
    mod.set_uplift( this_int_map["uplift_mode"], this_float_map["spinup_U"] );
    mod.set_current_frame(1);
    mod.run_components_combined();

    if(this_bool_map["staged_spinup"])
    {
      current_end_time = this_float_map["spinup_time"]+2*this_float_map["spinup_time"];
      mod.set_endTime(current_end_time);
      mod.set_K(this_float_map["spinup_K"]*0.1);
      mod.set_uplift( this_int_map["uplift_mode"], this_float_map["spinup_U"]*0.1 );
      mod.run_components_combined();
    }
  }

  //============================================================================
  // This wraps a number of functions to create a diamond-square based initial condition
  // that is then dissected.
  //============================================================================
  if(this_bool_map["diamond_square_spinup"])
  {
    int this_frame;
    
    cout << "I am running a spinup that uses a diamond square initial surface and" << endl;
    cout <<"then dissects the landscape for you." << endl;
    mod.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"]);
    // we need to check the feature order to make sure it is not bigger than the largest dimension of the raster
    int smallest_dimension;
    if( this_int_map["NRows"] <= this_int_map["NCols"])
    {
      smallest_dimension = this_int_map["NRows"];
    }
    else
    {
      smallest_dimension = this_int_map["NCols"];
    }
    // now get the largest possible dimension
    float lps = floor( log2( float(smallest_dimension) ) );
    int largest_possible_scale = int(lps);
    if( this_int_map["diamond_square_feature_order"] > largest_possible_scale )
    {
      this_int_map["diamond_square_feature_order"] = largest_possible_scale;
    }
    cout << "    The largest dimnesion is: " << smallest_dimension << " pixels." << endl;
    cout << "    So the biggest scale is: " << pow(2,largest_possible_scale) << " pixels or 2 to the " << largest_possible_scale << " power" << endl;
    cout << "    I am setting the scale to: " << this_int_map["diamond_square_feature_order"] << endl;

    // Now actually build the fractal surface
    mod.intialise_diamond_square_fractal_surface(this_int_map["diamond_square_feature_order"], this_float_map["diamond_square_relief"]);

    cout << "   Let me taper the edges of this fractal surface for you so that everything drains to the edge." << endl;
    cout << "   I am tapering along the " << this_int_map["taper_rows"] << " row closes to the N and S boundaries" << endl;
    mod.initialise_taper_edges_and_raise_raster(this_int_map["taper_rows"]);
    
    cout << "   I am superimposing a parabola with a relief of "  << this_float_map["parabola_relief"] << " metres" << endl;
    mod.superimpose_parabolic_surface(this_float_map["parabola_relief"]);
    
    cout << "   I am going to fill and then roughen the surface for you. Roughness elements will " << endl;
    cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
    mod.raise_and_fill_raster();
    mod.set_noise(this_float_map["roughness_relief"]);
    mod.random_surface_noise();
    mod.raise_and_fill_raster(); 
    
    this_frame = 9999;
    mod.print_rasters_and_csv( this_frame );
    
    
    cout << "Now I am going to run some fluvial incision at a small timestep" << endl;
    cout <<"I am going to use a moderate K value and uplift, since I don't" << endl;
    cout << "want to overexcavate the surface" << endl;
    mod.set_hillslope(false);
    //mod.set_timeStep( 1 );
    int this_uplift_mode = 0;     // block uplift
    float this_K = 0.001;
    float this_U = 0.0005;
    mod.set_K(this_K);
    
    mod.set_uplift( this_uplift_mode, this_U );
    this_frame = 1;
    for (int cycle = 1; cycle<=10; cycle++)
    {
      cout << "Cycle " << cycle << " of 10" << endl;
      mod.set_timeStep( float(cycle) );
      for (int i = 0; i<21; i++)
      {
        
        if (i%10 == 0)
        cout << "Initial dissection; i = " << i << " of 20" << endl;
        mod.fluvial_incision_with_uplift();

      }
      
      mod.raise_and_fill_raster();
      //this_frame++; 
      //mod.print_rasters_and_csv( this_frame );
    }
    
    // now print the model result
    this_frame = 9998;
    mod.print_rasters_and_csv( this_frame );
    
    current_end_time = 0;
  }



  //============================================================================
  // This cycles between version with hillslopes and without to speed up the initial
  // condition.  
  //============================================================================
  if(this_bool_map["cyclic_spinup"])
  {
    int this_frame;
    
    cout << "I am going to try to spin up the model by cycling between hillslope diffusion on and off."  << endl;
    cout << "First I need to ensure the raster is raised, filled and roughened."  <<endl;
    mod.raise_and_fill_raster();
    mod.set_noise(this_float_map["roughness_relief"]);
    mod.random_surface_noise();
    mod.raise_and_fill_raster(); 
    mod.set_print_interval(this_int_map["print_interval"]);
    mod.set_hillslope(false);
    
    // run a few cycles to get a network and then fill
    mod.set_timeStep( 1 );
    for (int i = 0; i<100; i++)
    {
      if (i%10 == 0)
      cout << "Initial dissection; i = " << i+1 << " of 100" << endl;
      mod.fluvial_incision_with_uplift();
    }
    mod.set_timeStep( this_float_map["spinup_dt"] );
    mod.raise_and_fill_raster(); 
    this_frame = 9998;
    mod.print_rasters_and_csv( this_frame );
      
    mod.set_current_frame(1);
    for(int i =0; i< this_int_map["spinup_cycles"]; i++)
    {
      cout << "++CYCLE NUMBER: "  << i << "+++++" << endl;
      
      // first we do a little bit of fluvial action
      current_end_time = current_end_time+this_float_map["spinup_time"];

      mod.set_endTime(current_end_time);
      mod.set_K(this_float_map["spinup_K"]);
      mod.set_uplift( this_int_map["uplift_mode"], this_float_map["spinup_U"] );
      mod.run_components_combined();

      // now do a bit more fluvial only but at a different parameter values
      cout << "A bit more fluvial at a different uplift rate" << endl;
      current_end_time = current_end_time+this_float_map["spinup_time"];
      mod.set_endTime(current_end_time);

      // logic for changing the K or U for the cycles. 
      if( this_bool_map["cycle_K"])
      {
        mod.set_K(this_float_map["spinup_K"]*this_float_map["cycle_K_factor"]);
      }
      if (this_bool_map["cycle_U"])
      {
        mod.set_uplift( this_int_map["uplift_mode"], this_float_map["spinup_U"]*this_float_map["cycle_U_factor"] );
      }
      mod.run_components_combined();
    }
  }


  if(this_bool_map["force_dissect"])
  {
    cout << "I am going to try to dissect your landscape by setting n = 1, having an rapid uplift rate, " << endl;
    cout << " setting a very high K, and running for a while" << endl;
    mod.raise_and_fill_raster(); 
    mod.set_n(1.0);
    mod.set_K(0.001);
    mod.set_uplift( 0,0.001 );
    mod.set_timeStep(100);
    mod.set_hillslope(false);
    int dissect_steps = this_int_map["force_dissect_steps"];
    for(int i = 0; i< dissect_steps; i++ )
    {
      if (i%100 == 0)
      {
        cout << "dissecting, step: " << i << " of " << dissect_steps << endl;
      }
      mod.fluvial_incision_with_uplift();
    }
    
    // set the print intervals and other parameters
    mod.set_m(this_float_map["m"]);
    mod.set_n(this_float_map["n"]);
    mod.set_uplift_mode(this_int_map["uplift_mode"]);
    mod.set_print_hillshade(this_bool_map["write_hillshade"]);
    mod.set_timeStep( this_float_map["dt"] );
    mod.set_print_interval(this_int_map["print_interval"]);
    mod.raise_and_fill_raster(); 
    
    // now print the model result
    int this_frame = 9997;
    mod.print_rasters_and_csv( this_frame );
  }



  //============================================================================
  // Logic for snapping to steady state using equation 4a from Mudd et al 2014
  // Every pixel is considered to be a channel. If you want hillslopes
  // you will need to run the model after this with the hillslopes turned on
  //============================================================================
  if(this_bool_map["snap_to_steady"])
  {
    cout << "I am going to snap the landscape to steady state. " << endl;
    cout << "The way this works is that chi is calculated and then the steady state" << endl;
    cout << " elevation is calculated using equation 4a from Mudd et al 2014 JGR-ES" << endl;
    cout << "This method assumes only fluvial erosion. If you want hillslopes you need" << endl;
    cout << "To turn them on and then let the model run to steady state. " << endl;
    
    // first ensure that the raster is filled and there are no below sea level nodes
    mod.raise_and_fill_raster(); 
    
    // now snap!
    float new_K = mod.fluvial_snap_to_steady_state_tune_K_for_relief(this_float_map["snapped_to_steady_uplift"], this_float_map["snapped_to_steady_relief"]);
    cout << "Getting a steady solution for a landscape with relief of " << this_float_map["snapped_to_steady_relief"]
         << " metres and uplift of " << this_float_map["snapped_to_steady_uplift"]*1000 << " mm per year." << endl;
    cout << "The new K is: " << new_K << endl;
    mod.set_K(new_K);
    
    if(this_bool_map["print_snapped_to_steady_frame"])
    {
      cout << "I am printing the snapped surface, increasing the frame count by 1." << endl;
      int this_frame = mod.get_current_frame();
      mod.print_rasters_and_csv( this_frame );
      this_frame++;
      mod.set_current_frame(this_frame);
      
    }
  }




  //============================================================================
  // Logic for a rudimentary steady forcing of uplift
  //============================================================================
  if(this_bool_map["rudimentary_steady_forcing"])
  {
    cout << "I am going to run some very basic steady forcing." << endl;
    cout << "Let me check the boundary conditions!" << endl;
    mod.print_boundary_conditions_to_screen();
    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    else
    {
      cout << "This forcing includes hillslope diffusion." << endl;
      mod.set_hillslope(true);
    }
    cout << "Let me run some steady forcing for you. " << endl;
    cout << "Starting with a K of: " << this_float_map["rudimentary_steady_forcing_K"] << endl;
    current_end_time = current_end_time+this_float_map["rudimentary_steady_forcing_time"];
    mod.set_timeStep( this_float_map["dt"] );
    mod.set_endTime(current_end_time);
    mod.set_K(this_float_map["rudimentary_steady_forcing_K"]);
    mod.set_uplift( this_int_map["uplift_mode"], this_float_map["rudimentary_steady_forcing_uplift"] );
    mod.run_components_combined();

  }



  
  //============================================================================
  // Logic for cyclic forcing of uplift
  //============================================================================
  if(this_bool_map["run_cyclic_forcing"])
  {
    cout << "I am going to force your model through a series of cycles, varying" << endl;
    cout << " either K or U." << endl;
    //int this_frame;

    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    else
    {
      cout << "This forcing includes hillslope diffusion." << endl;
      mod.set_hillslope(true);
    }
    
    // get the K value for the desired relief
    float first_cycle_K;
    if(this_bool_map["set_fixed_relief"])
    {
      cout << "I am calculating a K value that will get a relief of " << this_float_map["fixed_relief"] << " metres" << endl;
      cout << " for an uplift rate of " << this_float_map["baseline_U_for_cyclic"]*1000 << " mm/yr" << endl; 
      first_cycle_K = mod.fluvial_calculate_K_for_steady_state_relief(this_float_map["baseline_U_for_cyclic"],this_float_map["fixed_relief"]);
      cout << "The K value is: " << first_cycle_K << endl;
    }
    else
    {
      first_cycle_K = this_float_map["baseline_K_for_cyclic"];
    }

    // We need to run the model for a few timesteps at a short dt and then fill
    // to make sure no baselevel nodes are created
    mod.set_timeStep( 1 );
    mod.set_K(first_cycle_K);
    for (int i = 0; i<100; i++)
    {
      if (i%10 == 0)
      cout << "Initial dissection; i = " << i+1 << " of 100" << endl;
      mod.fluvial_incision_with_uplift();
    }
    mod.set_timeStep( this_float_map["spinup_dt"] );
    mod.raise_and_fill_raster(); 


    // Let the user know what you are doing
    if( this_bool_map["cycle_K"])
    {
      cout << "I am going to vary the baseline K by a factor of " << this_float_map["cycle_K_factor"] << endl;
    }
    if (this_bool_map["cycle_U"])
    {
      cout << "I am going to vary the baseline U by a factor of " << this_float_map["cycle_U_factor"] << endl;
    }

    // now for the model run
    mod.set_print_interval(this_int_map["print_interval"]);
    current_end_time = 0;
    mod.set_timeStep( this_float_map["cyclic_dt"] );
    
    for(int i =0; i< this_int_map["cyclic_cycles"]; i++)
    {
      cout << "++CYCLE NUMBER: "  << i << "+++++" << endl;
      
      // first we do a little bit of fluvial action
      current_end_time = current_end_time+this_float_map["cyclic_forcing_time"];

      mod.set_endTime(current_end_time);
      mod.set_K(first_cycle_K);
      mod.set_uplift( this_int_map["uplift_mode"], this_float_map["baseline_U_for_cyclic"] );
      mod.run_components_combined();

      // now do a bit more fluvial only but at a different parameter values
      cout << "A bit more fluvial at a different uplift rate" << endl;
      current_end_time = current_end_time+this_float_map["cyclic_forcing_time"];
      mod.set_endTime(current_end_time);

      // logic for changing the K or U for the cycles. 
      if( this_bool_map["cycle_K"])
      {
        mod.set_K(first_cycle_K*this_float_map["cycle_K_factor"]);
      }
      if (this_bool_map["cycle_U"])
      {
        mod.set_uplift( this_int_map["uplift_mode"], this_float_map["baseline_U_for_cyclic"]*this_float_map["cycle_U_factor"] );
      }
      mod.run_components_combined();
    }
  }


  //============================================================================
  // Logic for a random forcing of uplift and/or K
  //============================================================================
  if(this_bool_map["run_random_forcing"])
  {
    string ran_uplift_fname = OUT_DIR+OUT_ID+"_randomU.csv";
    ofstream Uout;
    Uout.open(ran_uplift_fname.c_str());
    Uout << "time,end_time,uplift_rate_m_yr" << endl;
    
    
    cout << "I am going to force your model through a series of cycles, varying" << endl;
    cout << " either K or U." << endl;
    
    // get the seed for the random forcing. 
    long seed = time(NULL);
    
    
    float time_gap = this_float_map["maximum_time_for_random_cycle"]-this_float_map["minimum_time_for_random_cycle"];
    float U_gap = this_float_map["maximum_U_for_random_cycle"]-this_float_map["minimum_U_for_random_cycle"];
    
    
    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    else
    {
      cout << "This forcing includes hillslope diffusion." << endl;
      mod.set_hillslope(true);
    }
    



    if(this_bool_map["set_fixed_relief"])
    {
      // get the K value for the desired relief
      float first_cycle_K;
      cout << "I am calculating a K value that will get a relief of " << this_float_map["fixed_relief"] << " metres" << endl;
      cout << " for an uplift rate of " << this_float_map["minimum_U_for_random"]*1000 << " mm/yr" << endl; 
      first_cycle_K = mod.fluvial_calculate_K_for_steady_state_relief(this_float_map["minimum_U_for_random"],this_float_map["fixed_relief"]);
      cout << "The K value is: " << first_cycle_K << endl;
      mod.set_K(first_cycle_K);
    }
    else
    {
      cout << "I am using the background K value." << endl;
      mod.set_K(this_float_map["background_K"]);
    }

    // now for the model run
    mod.set_print_interval(this_int_map["print_interval"]);
    current_end_time = 0;
    mod.set_timeStep( this_float_map["random_dt"] );
    
    for(int i =0; i< this_int_map["random_cycles"]; i++)
    {
      cout << "++CYCLE NUMBER: "  << i << "+++++" << endl;
      
      // get the time of this cycle
      float this_time = time_gap*ran3(&seed)+this_float_map["minimum_time_for_random_cycle"];
      current_end_time = current_end_time+this_time;
      mod.set_endTime(current_end_time);
      
      
      
      // now for the uplift
      float this_U = U_gap*ran3(&seed)+this_float_map["minimum_U_for_random_cycle"];
      current_end_time = current_end_time+this_time;
      
      Uout << this_time <<"," << current_end_time << "," << this_U << endl;
      cout << "Time of: " << this_time << " with U of " << this_U*1000 << " mm/yr." << endl;
      
      mod.set_uplift( this_int_map["uplift_mode"], this_U );
      mod.run_components_combined();
    }
    Uout.close();
  }

  //============================================================================
  // Logic for simulations with spatially varying uplift or K
  // There are a nubmer of flags here. One day I will try to organise the
  // flags to they are more consistent but today is not that day.
  // Just trying to get this thing working for the m/n paper at the moment. 
  //============================================================================
  if(this_bool_map["spatially_varying_forcing"])
  {
    // start by raising and filling the model
    mod.raise_and_fill_raster(); 
    
    cout << "I am running a simulation with spatially varying forcing." << endl;
    LSDRaster this_K_raster;
    LSDRaster this_U_raster;
    
    float this_max_K;
    float this_min_K;
    
    // Calculate or set the K parameter depending on your choices about the simulation
    if(this_bool_map["calculate_K_from_relief"])
    {
      cout << "I am calculating a K value that will get a relief of " << this_float_map["fixed_relief"] << " metres" << endl;
      cout << " for an uplift rate of " << this_float_map["min_U_for_spatial_var"]*1000 << " mm/yr" << endl; 
      this_min_K = mod.fluvial_calculate_K_for_steady_state_relief(this_float_map["min_U_for_spatial_var"],this_float_map["fixed_relief"]);
      this_max_K = this_min_K*this_float_map["spatial_K_factor"];
      cout << "The maximum K is: " << this_max_K << " and the minimum K is: " << this_min_K << endl;
    }
    else
    {
      cout << "I am using the maximum and minimum K values you have given me." << endl;
      this_max_K = this_float_map["spatially_varying_max_K"];
      this_min_K = this_float_map["spatially_varying_min_K"];
      cout << "The maximum K is: " << this_max_K << " and the minimum K is: " << this_min_K << endl;
    }
     
     
    if(this_bool_map["spatially_varying_K"])
    {
      cout << "I am going to vary K." << endl;
      if(this_bool_map["load_K_raster"])
      {
        
        string header = DATA_DIR+DEM_ID+"_KRaster.hdr";
        cout << "The full read path is: " << header << endl;
        ifstream file_info_in;
        file_info_in.open(header.c_str());
        // check if the parameter file exists
        if( not file_info_in.fail() )
        {
          cout << "I found the header. I am loading this initial file. " << endl;
          LSDRaster temp_raster(DATA_DIR+DEM_ID+"_KRaster","bil");
          this_K_raster = temp_raster;
        }
        else
        {
          // If you can't read the file then turn the creation routine on
          this_bool_map["load_K_raster"] = false;
          cout << "Warning, the K raster you wanted to load doesn't exist. I am making a new one." << endl;
        }
      }
      if( not this_bool_map["load_K_raster"])
      {
        switch (this_int_map["spatial_K_method"])
        {
          case 0:
            {
              cout << "K variation case 0. Min K: " << this_min_K<<  " and max K: " << this_max_K << endl;
              LSDRasterMaker KRaster1(this_int_map["NRows"],this_int_map["NCols"]);
              KRaster1.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_min_K);
              KRaster1.random_square_blobs(this_int_map["min_blob_size"], this_int_map["max_blob_size"], 
                                    this_min_K, this_max_K,
                                    this_int_map["n_blobs"]);
                                    
              // smooth the raster
              //cout << "I am going to smooth K a few times" << endl;
              for(int si = 0; si< this_int_map["K_smoothing_steps"]; si++)
              {
                //cout << "diffuse_K step: " << si << endl;
                KRaster1.smooth(0);
              }
              this_K_raster = KRaster1.return_as_raster();
              
            }
            break;
          case 1:
            {
              cout << "HAHAHA This is a secret kill switch. You lose! Try again next time Sonic!"  << endl;
              exit(EXIT_FAILURE);
            }
            break;
          default:
            {
              cout << "K variation. The options are 0 == random squares" << endl;
              cout << "  0 == random squares" << endl;
              cout << "  1 == sine waves (I lied, at the moment this doesn't work--SMM Sept 2017)." << endl;
              cout << "You didn't choose a valid option so I am defaulting to random squares." << endl;
              cout << "K variation: Min K: " << this_min_K<<  " and max K: " << this_max_K << endl;
              LSDRasterMaker KRaster1(this_int_map["NRows"],this_int_map["NCols"]);
              KRaster1.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_min_K);
              KRaster1.random_square_blobs(this_int_map["min_blob_size"], this_int_map["max_blob_size"], 
                                    this_min_K, this_max_K,
                                    this_int_map["n_blobs"]);
  
              // smooth the raster
              for(int si = 0; si< this_int_map["K_smoothing_steps"]; si++)
              {
                //cout << "diffuse_K step: " << si << endl;
                KRaster1.smooth(0);
              }
              this_K_raster = KRaster1.return_as_raster();
            }
            break;
        }
      }
    }
    else
    {
      cout<< "You decided not to vary K. The value of K is: " << this_min_K << endl;
      LSDRasterMaker KRaster2(this_int_map["NRows"],this_int_map["NCols"]);
      KRaster2.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_min_K);
      this_K_raster = KRaster2.return_as_raster();
    }
    
    
    
    // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_KRaster";
    string bil_name = "bil";
    this_K_raster.write_raster(K_fname,bil_name);
    
    // now deal with the uplift. We only use a sine uplift for now. 
    if(this_bool_map["spatially_varying_U"])
    {
      cout << "I am varying U." << endl;
      if(this_bool_map["load_U_raster"])
      {
        
        string header = DATA_DIR+DEM_ID+"_URaster.hdr";
        cout << "The full read path is: " << header << endl;
        ifstream file_info_in;
        file_info_in.open(header.c_str());
        // check if the parameter file exists
        if( not file_info_in.fail() )
        {
          cout << "I found the header. I am loading this initial file. " << endl;
          LSDRaster temp_raster(DATA_DIR+DEM_ID,"bil");
          this_U_raster = temp_raster;
        }
        else
        {
          // If you can't read the file then turn the creation routine on
          this_bool_map["load_U_raster"] = false;
          cout << "Warning, the U raster you wanted to load doesn't exist. I am making a new one." << endl;
        }
      }
      if( not this_bool_map["load_U_raster"])
      {
        switch (this_int_map["spatial_U_method"])
        {
          case 0:
            {
              cout << "Spatial variation in U. HAHAHA This is a secret kill switch. You lose! Try again next time Sonic!"  << endl;
              exit(EXIT_FAILURE);
            }
          case 1:
            {
              cout << "Spatial variation in U. Case 0." << endl;
              LSDRasterMaker URaster(this_int_map["NRows"],this_int_map["NCols"]);
              URaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["min_U_for_spatial_var"]);
              
              // The way the sine function work is that you give the function 
              // coefficients that give the varuous amplitudes of sine waves with 
              // 1/2 wavelengths that are 1*model_domain, 1/2*model_domain, 1/3*model domain, etc.
              // Here we only do the y coefficients since we are only going to vary
              // uplift in the y direction.
              vector<float> x_coeff;
              vector<float> y_coeff;
              y_coeff.push_back(10);
              URaster.sine_waves(x_coeff, y_coeff);
              
              // now scale to the desired minimum and maximum
              URaster.scale_to_new_minimum_and_maximum_value(float_default_map["min_U_for_spatial_var"],
                                                             float_default_map["max_U_for_spatial_var"]);
              this_U_raster = URaster.return_as_raster();
            }
            break;
  
  
            break;
          default:
            {
              cout << "Spatial variation in U. The options are 0 == random squares" << endl;
              cout << "  0 == random squares (I lied, at the moment this doesn't work--SMM Sept 2017)." << endl;
              cout << "  1 == sine waves" << endl;
              cout << "You didn't choose a valid option so I am defaulting to random squares." << endl;
              LSDRasterMaker URaster(this_int_map["NRows"],this_int_map["NCols"]);
              URaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["min_U_for_spatial_var"]);
              
              // The way the sine function work is that you give the function 
              // coefficients that give the varuous amplitudes of sine waves with 
              // 1/2 wavelengths that are 1*model_domain, 1/2*model_domain, 1/3*model domain, etc.
              // Here we only do the y coefficients since we are only going to vary
              // uplift in the y direction.
              vector<float> x_coeff;
              vector<float> y_coeff;
              y_coeff.push_back(10);
              URaster.sine_waves(x_coeff, y_coeff);
              
              // now scale to the desired minimum and maximum
              URaster.scale_to_new_minimum_and_maximum_value(float_default_map["min_U_for_spatial_var"],
                                                             float_default_map["max_U_for_spatial_var"]);
              this_U_raster = URaster.return_as_raster();
            }
            break;
        }
      }
    }
    else
    {
      cout << "I am not going to vary U in space. Instead I will use block uplift of " << this_float_map["min_U_for_spatial_var"] << endl;
      LSDRasterMaker URaster(this_int_map["NRows"],this_int_map["NCols"]);
      URaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["min_U_for_spatial_var"]);
      this_U_raster = URaster.return_as_raster();
    }
    // write the raster
    string U_fname = OUT_DIR+OUT_ID+"_URaster";
    this_U_raster.write_raster(U_fname,bil_name);
    
    
    cout << "========================================" << endl;
    cout << "I am now going to run your model with some spatially variable parameters." << endl;
    
    // set the end time and other model parameters
    if( not this_bool_map["hillslopes_on"] )
    {
      cout << "I'm turning hillslope diffusion off." << endl;
      mod.set_hillslope(false);
    }
    else
    {
      cout << "This forcing includes hillslope diffusion." << endl;
      mod.set_hillslope(true);
    }
    
    bool use_adaptive_timestep =   this_bool_map["use_adaptive_timestep"];
    
    if (use_adaptive_timestep)
    {
      cout << "I am using an adaptive timestep" << endl;
    }
    else
    {
      cout << "I am using a fixed timestep" << endl;
    }
    
    mod.set_maxtimeStep(this_float_map["maximum_timestep"]);
    mod.set_print_interval(this_int_map["print_interval"]);
    mod.set_timeStep( this_float_map["spatial_dt"] );
    
    mod.set_next_printing_time(0);
    mod.set_float_print_interval(this_float_map["float_print_interval"]);
    
    if(this_bool_map["snap_to_steep_for_spatial_uplift"])
    {
      
      cout << "I am going to snap this model to a steep landscape. " << endl;
      mod.set_K(this_K_raster.get_data_element(0,0));
      
      if(this_bool_map["snap_to_minimum_uplift"])
      {
        cout << "I'm snapping to the minimum uplift rate." << endl;
        mod.fluvial_snap_to_steady_state(this_float_map["min_U_for_spatial_var"]);
      }
      else
      {
        cout << "I'm snapping to the maximum uplift rate." << endl;
        mod.fluvial_snap_to_steady_state(this_float_map["max_U_for_spatial_var"]);
      }
      int snap_frame = 9991;
      mod.print_rasters_and_csv( snap_frame );
    }
    
    // Now run the model
    // Use cycles and fill after to avoid internal baselevel nodes
    for(int i = 0; i< this_int_map["spatial_cycles"]; i++)
    {
      current_end_time = current_end_time+this_float_map["spatial_variation_time"];
      mod.set_endTime(current_end_time);
      
      mod.run_components_combined(this_U_raster, this_K_raster,use_adaptive_timestep);
      mod.raise_and_fill_raster(); 
    }
    
    // now run at steady condition for a few extra cycles
    current_end_time = current_end_time+float(this_int_map["spatial_cycles"])*this_float_map["spatial_variation_time"];
    mod.set_endTime(current_end_time);
    mod.run_components_combined(this_U_raster, this_K_raster,use_adaptive_timestep);
  }
  
}