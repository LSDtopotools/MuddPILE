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
#include "../LSDSpatialCSVReader.hpp"
#include "../LSDLithoCube.hpp"
#include "../LSDRaster.hpp"
#include "../LSDIndexRaster.hpp"
using namespace std;


int main (int nNumberofArgs,char *argv[])
{

  cout << "================================================================" << endl;
  cout << "|| Welcome to the MuddPILEdirver!                             ||" << endl;
  cout << "|| This program drives the MuddPILE lanscape evolution model. ||" << endl;
  cout << "|| One day this model will have documentation.                ||" << endl;
  cout << "|| This program was developed by Simon M. Mudd                ||" << endl;
  cout << "|| James Jenkinson, Declan Valters and Fiona Clubb            ||" << endl;
  cout << "||  at the University of Edinburgh                            ||" << endl;
  cout << "================================================================" << endl;
  cout << "|| If you use these routines please cite:                     ||" << endl;
  cout << "|| https://www.doi.org/10.5281/zenodo.997388                  ||" << endl;
  cout << "|| The initial paper using this model is:                     ||" << endl;
  cout << "|| https://onlinelibrary.wiley.com/doi/full/10.1002/esp.3923  ||" << endl;
  cout << "================================================================" << endl;
  cout << "|| Documentation can be found at:                             ||" << endl;
  cout << "|| https://lsdtopotools.github.io/LSDTT_documentation/        ||" << endl;
  cout << "================================================================" << endl;

  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);


  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

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

  // Some logi for filling rasters that are being loaded
  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  bool_default_map["raster_is_filled"] = false;
  bool_default_map["print_fill_raster"] = false;



  // Read a channel file
  bool_default_map["extract_single_channel"] = false;
  string_default_map["channel_source_fname"] = "single_channel_source";
  string_default_map["single_channel_print_fname"] = "single_channel";

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

  // More complex snapping
  bool_default_map["snap_to_steady_variable_variable_K_variable_U"] = false;
  string_default_map["variable_K_name"] = "variable_K";
  string_default_map["variable_U_name"] = "variable_U";
  string_default_map["fixed_channel_csv_name"] = "single_channel_nodes";  // This contains the fixed channel data

  // Even more complex snapping
  bool_default_map["snap_to_steady_variable_variable_K_variable_U_use_LithoCube"] = false;

  // Yet more complex snapping that attempts to migrate divides by cycling snapping
  bool_default_map["cycle_fluvial_snapping"] = false;
  int_default_map["snapping_cycles"] = 30;
  int_default_map["snap_print_rate"] = 10;     // how frequently the snapping is printed in terms of snap cycles

  // Yet more complex snapping that attempts to migrate divides by cycling snapping of the critical slope component
  bool_default_map["cycle_hillslope_snapping"] = false;
  int_default_map["snapping_hillslope_cycles"] = 1;
  int_default_map["snap_hillslope_print_rate"] = 10;     // how frequently the snapping is printed in terms of snap cycles  

  // This is for the critical slope
  float_default_map["rudimentary_critical_slope"] = 0.8;
  int_default_map["threshold_contributing_pixels"] = 1000;
  bool_default_map["snap_to_critical_slopes"] = false;
  string_default_map["Sc_raster_fname"] = "S_c";
  bool_default_map["use_Sc_raster"] = false;
  bool_default_map["snap_to_steady_critical_slopes_use_LithoCube"] = false;
  bool_default_map["print_Sc_raster"] = false; 

  // Smoothing parameters 
  float_default_map["smoothing_window_radius"] =50;
  int_default_map["smoothing_sweeps"] = 1;

  // forcing of dissection
  bool_default_map["force_dissect"] = true;
  int_default_map["force_dissect_steps"] = 10000;

  // Some parameters for very rudimentary steady forcing
  bool_default_map["make_constant_K_raster"] = false;
  bool_default_map["make_constant_U_raster"] = false;
  bool_default_map["rudimentary_steady_forcing"] = false;
  float_default_map["rudimentary_steady_forcing_time"] = 100000;
  float_default_map["rudimentary_steady_forcing_uplift"] = 0.0005;
  float_default_map["rudimentary_steady_forcing_K"] = 0.0005;

  // Parameters for hillslopes
  bool_default_map["hillslopes_on"] = false;
  bool_default_map["nonlinear"] = false;

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
  bool_default_map["make_spatially_varying_K"] = false;   // This just prints a raster
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
  bool_default_map["make_spatially_varying_U"] = false;  // This just prints a raster
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
  int_default_map["increase_amt"] = 10;
  bool_default_map["finish_with_steady_forcing"] = false;

  bool_default_map["simple_fault_uplift"] = false;
  int_default_map["uplift_amt"] = 10;
  int_default_map["uplift_time"] = 50000;
  float_default_map["uplift_rate"] = 0.0002;
  bool_default_map["base_level_fall"] = false;
  int_default_map["transient_print_interval"] = 2;

  // tilt scenarios
  bool_default_map["run_instantaneous_tilt"] = false;
  bool_default_map["run_progressive_tilt"] = false;
  int_default_map["tilt_time"] = 100000;
  float_default_map["tilt_angle"] = 10;
  string_default_map["tilt_boundary"] = "S";
  int_default_map["tilt_print_interval"] = 150;

  // Various options to load and process lithologic data
  bool_default_map["load_lithology"] = false;
  bool_default_map["only_test_lithology"] = false;
  string_default_map["vo_filename"] = "test.vo";
  string_default_map["vo_asc_filename"] = "test.asc";
  string_default_map["lookup_filename"] = "lookup_table_strati_K.csv";
  bool_default_map["print_K_raster"] = false;
  bool_default_map["print_lithocode_raster"] = false;


  //Option to adjust elevation of raster before lithocube snapping
  int_default_map["elevation_change"] = 0;


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
  // Logic for loading lithology
  //============================================================================ 
  if(this_bool_map["only_test_lithology"])
  {
    cout << "Let's load the lithology!";
    string this_vo_filename = DATA_DIR+this_string_map["vo_filename"];
    cout << "I am going to load lithology information" << endl;
    cout << "The filename is: " << this_vo_filename << endl;
    LSDLithoCube LSDLC(this_vo_filename);

    string this_vo_asc = DATA_DIR+this_string_map["vo_asc_filename"];
    LSDLC.ingest_litho_data(this_vo_asc);

    // Print a layer to a raster to see if everything is okay
    // I am going to need to impose the georeferencing


    LSDRaster elev_raster(DATA_DIR+DEM_ID,"bil");
    
    // Deal with the georeferencing
    int UTM_zone;
    bool is_North = true;
    string NorS = "N";
    elev_raster.get_UTM_information(UTM_zone,is_North);
    if (is_North == false)
    {
      NorS = "S";
    }  
    LSDLC.impose_georeferencing_UTM(UTM_zone, NorS);

    

    // test ingesting strati to K lookup table
    // cout << "I will now ingest the lookup table.";
    // string this_lookup_table = DATA_DIR+this_string_map["lookup_filename"];
    // LSDLC.read_strati_to_K_csv(this_lookup_table);

    LSDLC.hardcode_fill_strati_to_K_map();

    // // test printing of a layer
    // string test_layer_name = DATA_DIR+DEM_ID+"_Layer61";
    // int test_layer = 61;
    // LSDIndexRaster test_layer_raster = LSDLC.get_raster_from_layer(test_layer);
    // LSDRaster K_raster_from_lithocube_layer = LSDLC.index_to_K(test_layer_raster, 0.000001);
    // cout << "Printing a K raster from layer " << test_layer << endl;
    // K_raster_from_lithocube_layer.write_raster(test_layer_name, "bil");

    string lithocodes_raster_name = DATA_DIR+OUT_ID+"_LithoCodes";
    string test_K_raster_name = DATA_DIR+OUT_ID+"_KRasterLSDLC";

    // LSDLC.write_litho_raster_from_elevations(elev_raster, test_litho_raster_name, "bil");
    cout << "Getting the lithology codes" << endl;
    LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster);
    LSDRaster lithocodes_raster(lithocodes_index_raster);
    cout << "Printing a lithocodes raster" << endl;
    lithocodes_raster.write_raster(lithocodes_raster_name, "bil");    
    cout << "Converting the litho codes to K values" << endl;
    LSDRaster K_raster_from_lithocube = LSDLC.index_to_K(lithocodes_index_raster, 0.000001);
    cout << "Printing a K raster" << endl;
    K_raster_from_lithocube.write_raster(test_K_raster_name, "bil");
    cout << "Wrote a K raster for you" << endl;

    //LSDLC.get_K_raster_from_raster(elev_raster);



    cout << "I've tested the lithology object and now am exiting." << endl;
    exit(0);
    

  } 

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
    cout << "The full read path for the header (so I can see if the raster exists) is: " << header << endl;
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

      map<string,string> GRS = temp_raster.get_GeoReferencingStrings();
      string CS = GRS["ENVI_coordinate_system"];
      string MI = GRS["ENVI_map_info"];
      cout << "I am laoding a raster with Cooridinate string: " << endl;
      cout << CS << endl;
      cout << "and MI is: " << endl;
      cout << MI << endl;

      this_int_map["NRows"] = temp_raster.get_NRows();
      this_int_map["NCols"] = temp_raster.get_NCols();
      this_float_map["DataResolution"] = temp_raster.get_DataResolution();


      //==========================================================================
      // Fill the raster
      //==========================================================================
      LSDRaster filled_topography,carved_topography;
      if ( this_bool_map["raster_is_filled"] )
      {
        cout << "You have chosen to use a filled raster." << endl;
        filled_topography = temp_raster;
      }
      else
      {
        cout << "Let me fill that raster for you, the min slope is: "
            << this_float_map["min_slope_for_fill"] << endl;
        if(this_bool_map["carve_before_fill"])
        {
          carved_topography = temp_raster.Breaching_Lindsay2016();
          filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
        }
        else
        {
          filled_topography = temp_raster.fill(this_float_map["min_slope_for_fill"]);
        }
      }

      // Print the fill raster if you want it
      if (this_bool_map["print_fill_raster"])
      {
        cout << "Let me print the fill raster for you."  << endl;
        string filled_raster_name = OUT_DIR+OUT_ID+"_Fill";
        filled_topography.write_raster(filled_raster_name,raster_ext);
      }

      // This gets a single channel from your loaded raster for use later
      if(this_bool_map["extract_single_channel"])
      {

        LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

        LSDRasterInfo RI(filled_topography);
        // Get the latitude and longitude
        cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["channel_source_fname"]) );

        // Get the local coordinates
        vector<float> fUTM_easting,fUTM_northing;
        source_points_data.get_x_and_y_from_latlong(fUTM_easting,fUTM_northing);
        float X,Y;

        // This grabs the single channel
        if ( int(fUTM_easting.size()) != 0)
        {

          X = fUTM_easting[0];
          Y = fUTM_northing[0];

          cout << "Looking for the single source. " << endl;
          cout << "Easting is is: " << X << " and northing is: " << Y << endl;


          vector<int> node_list = FlowInfo.get_flow_path(X,Y);

          cout << "Your node list has " << node_list.size() << " nodes." << endl;

          LSDRaster FD = FlowInfo.distance_from_outlet();
          LSDRaster DA = FlowInfo.write_DrainageArea_to_LSDRaster();

          string fname = this_string_map["single_channel_print_fname"];
          cout << "Printing your channel to a csv file." << endl;
          FlowInfo.print_vector_of_nodeindices_to_csv_file_with_latlong(node_list,DATA_DIR, fname,
                                                        filled_topography, FD, DA);


        }
        else
        {
          cout << "You are trying to make a channel but I cannot find the source." << endl;
        }

      }

      LSDRasterModel temp_mod(temp_raster);
      mod = temp_mod;
      // taper the edges of the input raster to 0 at N and S boundaries
      mod.initialise_taper_edges_and_raise_raster(1);
      create_initial_surface = false;
      this_bool_map["force_dissect"] = false;
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
  else
  {
    cout << "I am not creating an initial surface. This is probably because you chose to load a raster. " << endl;
  }


  //============================================================================
  // Logic for making spatially varying K fields
  //============================================================================
  if(this_bool_map["make_spatially_varying_K"])
  {
    LSDRaster Temp_raster = mod.return_as_raster();
    cout << "I am going to make a spatially varying raster for the K parameter." << endl;

    LSDRasterMaker KRaster(Temp_raster);
    KRaster.set_to_constant_value(this_float_map["spatially_varying_min_K"]);
    KRaster.random_square_blobs(this_int_map["min_blob_size"], this_int_map["max_blob_size"],
                                this_float_map["spatially_varying_min_K"], this_float_map["spatially_varying_max_K"],
                                this_int_map["n_blobs"]);

    // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_KRaster";
    string bil_name = "bil";

    KRaster.write_raster(K_fname,bil_name);
  }

  if(this_bool_map["make_constant_K_raster"])
  {
    LSDRaster Temp_raster = mod.return_as_raster();
    cout << "I am going to make a contant raster for the K parameter." << endl;
    cout << "Setting K to "  << this_float_map["rudimentary_steady_forcing_K"] << endl;
    LSDRasterMaker KRaster(Temp_raster);
    KRaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_K"]);

    // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_ConstKRaster";
    string bil_name = "bil";

    KRaster.write_raster(K_fname,bil_name);
  }

  if(this_bool_map["make_spatially_varying_U"])
  {
    // now do a sine version.
    LSDRaster Temp_raster = mod.return_as_raster();
    cout << "I am going to make a spatially varying U raster. This just prints the raster! It only uses the sine version." << endl;
    LSDRasterMaker URaster(Temp_raster);
    URaster.set_to_constant_value(this_float_map["min_U_for_spatial_var"]);
    vector<float> x_coeff;
    vector<float> y_coeff;
    y_coeff.push_back(10);
    //y_coeff.push_back(0);
    //y_coeff.push_back(2);
    //y_coeff.push_back(0);
    //y_coeff.push_back(5);

    URaster.sine_waves(x_coeff, y_coeff);

    // now scale to the desired minimum and maximum
    URaster.scale_to_new_minimum_and_maximum_value(float_default_map["min_U_for_spatial_var"],
                                                    float_default_map["max_U_for_spatial_var"]);
    string U_fname = OUT_DIR+OUT_ID+"_URaster";
    string bil_name = "bil";
    URaster.write_raster(U_fname,bil_name);
  }

  if(this_bool_map["make_constant_U_raster"])
  {
    LSDRaster Temp_raster = mod.return_as_raster();
    cout << "I am going to make a contant raster for the U parameter." << endl;

    LSDRasterMaker KRaster(Temp_raster);
    KRaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_uplift"]);

    // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_ConstURaster";
    string bil_name = "bil";

    KRaster.write_raster(K_fname,bil_name);
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
  else
  {
    cout << "You will not spin up the surface." << endl;
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
  else
  {
    cout << "We are not doing the diamond square spinup routines." << endl;
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

  if(this_bool_map["snap_to_steady_variable_variable_K_variable_U"])
  {
    cout << "I am going to snap to steady using various variable K and U fields." << endl;
    cout << "I will also use a fixed channel if you have one. " << endl;

    // First load the two variable K and variable U maps
    string U_fname = DATA_DIR+this_string_map["variable_U_name"];
    string U_header = DATA_DIR+this_string_map["variable_U_name"]+".hdr";
    string K_fname = DATA_DIR+this_string_map["variable_K_name"];
    string K_header = DATA_DIR+this_string_map["variable_K_name"]+".hdr";

    LSDRaster Temp_raster = mod.return_as_raster();

    LSDRaster K_values;
    ifstream K_file_info_in;
    K_file_info_in.open(K_header.c_str());
    // check if the parameter file exists
    if( not K_file_info_in.fail() )
    {
      LSDRaster temp_K_values(K_fname,raster_ext);
      K_values = temp_K_values;
    }
    else
    {
      cout << "No variable K raster found. I will make one with a constant value of: " << this_float_map["rudimentary_steady_forcing_K"] << endl;
      LSDRasterMaker KRaster(Temp_raster);
      cout << "Let me set the value." << endl;
      KRaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_K"]);
      cout << "Set to constant K values. " << endl;
      K_values = KRaster.return_as_raster();
    }
    K_file_info_in.close();

    LSDRaster U_values;
    ifstream U_file_info_in;
    U_file_info_in.open(U_header.c_str());
    // check if the parameter file exists
    if( not U_file_info_in.fail() )
    {
      LSDRaster temp_U_values(U_fname,raster_ext);
      U_values = temp_U_values;
    }
    else
    {
      cout << "No variable U raster found. I will make one with a constant value of " << this_float_map["rudimentary_steady_forcing_uplift"] << endl;
      LSDRasterMaker URaster(Temp_raster);
      URaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_uplift"]);
      U_values = URaster.return_as_raster();
    }
    U_file_info_in.close();


    LSDRasterInfo RI(K_values);
    // Get the latitude and longitude
    cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
    LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

    int n_snaps;
    if (not this_bool_map["cycle_fluvial_snapping"])
    {
      n_snaps = 1;
    }
    else
    {
      n_snaps = this_int_map["snapping_cycles"];
    }

    int this_frame = 5000;
    mod.set_current_frame(this_frame);
    cout << "Let me do a few snaps. The number of snaps is: " << n_snaps << endl;
    for (int i = 0; i < n_snaps; i++)
    {
      // Now the actual snapping!
      cout << endl << endl << endl << endl;
      cout << "=================================================" << endl;
      cout << "The U value at (500,500) is" << U_values.get_data_element(500,500) << endl;
      mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"]);

      // Now print

      if (i % this_int_map["snap_print_rate"] == 0)
      {
        cout << "This is snap " << i+1 << " of " << n_snaps << endl;

        if(this_bool_map["print_snapped_to_steady_frame"])
        {
          cout << endl << endl << "====================================" << endl;
          cout << "I am printing the snapped surface, increasing the frame count by 1." << endl;
          cout << "The frame number is " << this_frame << endl;
          mod.set_current_frame(this_frame);
          mod.print_rasters_and_csv( this_frame );
        }
      }
      this_frame++;

    }

    cout << "Thanks for snapping!" << endl;
  }























if(this_bool_map["snap_to_steady_variable_variable_K_variable_U_use_LithoCube"])
  {
    cout << "I am going to snap to steady using various variable K and U fields, and 3D lithological data." << endl;
    cout << "I will also use a fixed channel if you have one. " << endl;

    // First load the two variable K and variable U maps
    string U_fname = DATA_DIR+this_string_map["variable_U_name"];
    string U_header = DATA_DIR+this_string_map["variable_U_name"]+".hdr";

    // Then load the vo and asc files containing the 3D litho data 
    string this_vo_filename = DATA_DIR+this_string_map["vo_filename"];
    string this_vo_asc = DATA_DIR+this_string_map["vo_asc_filename"];

    // Set names for K raster and lithocode raster to be printed 
    string K_name = DATA_DIR+OUT_ID+"_KRasterLSDLC";
    string lithocodes_raster_name = DATA_DIR+OUT_ID+"_LithoCodes";

    // Load the elevation raster
    LSDRaster elev_raster(DATA_DIR+DEM_ID,"bil");
    
    // Adjust the elevation of the raster
    if (this_int_map["elevation_change"] != 0)
    {
      elev_raster.AdjustElevation(this_int_map["elevation_change"]);
      cout << "Adjusted the elevation of your raster by " << this_int_map["elevation_change"] << endl;
      elev_raster.write_raster(DATA_DIR+DEM_ID+"_Dropped", "bil");
      // I will also want to adjust the elevation of the fixed channel here so I don't have to edit the csv file manually
    }
    
    // Initialise the lithocube and ingest the lithology data
    cout << "I am going to load lithology information" << endl;
    cout << "The filename is: " << this_vo_filename << endl;
    LSDLithoCube LSDLC(this_vo_filename);
    LSDLC.ingest_litho_data(this_vo_asc);    
    
    // Deal with the georeferencing
    int UTM_zone;
    bool is_North = true;
    string NorS = "N";
    elev_raster.get_UTM_information(UTM_zone,is_North);
    if (is_North == false)
    {
      NorS = "S";
    }  
    LSDLC.impose_georeferencing_UTM(UTM_zone, NorS);

    // Let's create the map of stratigraphy codes and K values 
    LSDLC.hardcode_fill_strati_to_K_map();    

    // Now we'll get the lithology codes, convert them to K values and make a raster of K values for the model to use
    cout << "Getting the lithology codes" << endl;
    LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster);
    cout << "Converting the litho codes to K values" << endl;
    LSDRaster K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000001);
       
    LSDRaster Temp_raster = mod.return_as_raster();

    // Make a raster of uplift values 
    LSDRaster U_values;
    ifstream U_file_info_in;
    U_file_info_in.open(U_header.c_str());
    // check if the parameter file exists
    if( not U_file_info_in.fail() )
    {
      LSDRaster temp_U_values(U_fname,raster_ext);
      U_values = temp_U_values;
    }
    else
    {
      cout << "No variable U raster found. I will make one with a constant value of " << this_float_map["rudimentary_steady_forcing_uplift"] << endl;
      LSDRasterMaker URaster(Temp_raster);
      URaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_uplift"]);
      U_values = URaster.return_as_raster();
    }
    U_file_info_in.close();


    LSDRasterInfo RI(K_values);
    // Get the latitude and longitude
    cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
    LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

    int n_snaps;
    if (not this_bool_map["cycle_fluvial_snapping"])
    {
      n_snaps = 1;
    }
    else
    {
      n_snaps = this_int_map["snapping_cycles"];
    }

    int this_frame = 5000;
    mod.set_current_frame(this_frame);
    cout << "Let me do a few snaps. The number of snaps is: " << n_snaps << endl;
    for (int i = 0; i < n_snaps; i++)
    {
      // Now the actual snapping!
      cout << endl << endl << endl << endl;
      cout << "=================================================" << endl;
      cout << "The U value at (500,500) is" << U_values.get_data_element(500,500) << endl;

      cout << "*****************************************" << endl;
      cout << "This is snap " << i+1 << " of " << n_snaps << endl;
      cout << "*****************************************" << endl;
      
      LSDRaster lithocodes_raster;

      // if(i>0)
      // {
      cout << "Let me create a raster of the topography from the previous loop." << endl;
      LSDRaster new_muddpile_topo_raster = mod.return_as_raster();
      cout << "Now I'll check where in the lithocube we are and make a raster of lithocodes." << endl;
      LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(new_muddpile_topo_raster);
      cout << "I am now converting the lithocodes index raster to an LSDRaster" << endl;
      lithocodes_raster = LSDRaster(lithocodes_index_raster);
      cout << "Converting the litho codes to K values" << endl;
      K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000001);

      // }

      mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"]);

      // Now print

      if (i % this_int_map["snap_print_rate"] == 0)
      {
        cout << "This is snap " << i+1 << " of " << n_snaps << endl;

        if(this_bool_map["print_snapped_to_steady_frame"])
        {
          cout << endl << endl << "====================================" << endl;
          cout << "I am printing the snapped surface, increasing the frame count by 1." << endl;
          cout << "The frame number is " << this_frame << endl;
          mod.set_current_frame(this_frame);
          mod.print_rasters_and_csv( this_frame );
          if(this_bool_map["print_K_raster"])
          {
            cout << "Printing a K raster" << endl;
            K_values.write_raster(K_name + to_string(this_frame), "bil");
          }
          if(this_bool_map["print_lithocode_raster"]) //&& i>0
          {
            cout << "Printing a lithocode raster" << endl;
            lithocodes_raster.write_raster(lithocodes_raster_name + to_string(this_frame), "bil");
          }

        }
      }
      this_frame++;

    }

    cout << "Thanks for snapping!" << endl;
  }









if(this_bool_map["snap_to_steady_critical_slopes_use_LithoCube"])
{
  cout << "I am going to snap to steady using various variable K and U fields, and 3D lithological data." << endl;
  cout << "I will also use a fixed channel if you have one. " << endl;

  // First load the two variable K and variable U maps
  string U_fname = DATA_DIR+this_string_map["variable_U_name"];
  string U_header = DATA_DIR+this_string_map["variable_U_name"]+".hdr";

  // Then load the vo and asc files containing the 3D litho data 
  string this_vo_filename = DATA_DIR+this_string_map["vo_filename"];
  string this_vo_asc = DATA_DIR+this_string_map["vo_asc_filename"];

  // Set names for K raster and lithocode raster to be printed 
  string K_name = DATA_DIR+OUT_ID+"_KRasterLSDLC";
  string Sc_name = DATA_DIR+OUT_ID+"_ScRasterLSDLC";
  string lithocodes_raster_name = DATA_DIR+OUT_ID+"_LithoCodes";

  // Load the elevation raster
  LSDRaster elev_raster(DATA_DIR+DEM_ID,"bil");
  
  // Adjust the elevation of the raster
  if (this_int_map["elevation_change"] != 0)
  {
    elev_raster.AdjustElevation(this_int_map["elevation_change"]);
    cout << "Adjusted the elevation of your raster by " << this_int_map["elevation_change"] << endl;
    elev_raster.write_raster(DATA_DIR+DEM_ID+"_Dropped", "bil");
    // I will also want to adjust the elevation of the fixed channel here so I don't have to edit the csv file manually
  }
  
  // Initialise the lithocube and ingest the lithology data
  cout << "I am going to load lithology information" << endl;
  cout << "The filename is: " << this_vo_filename << endl;
  LSDLithoCube LSDLC(this_vo_filename);
  LSDLC.ingest_litho_data(this_vo_asc);    
  
  // Deal with the georeferencing
  int UTM_zone;
  bool is_North = true;
  string NorS = "N";
  elev_raster.get_UTM_information(UTM_zone,is_North);
  if (is_North == false)
  {
    NorS = "S";
  }  
  LSDLC.impose_georeferencing_UTM(UTM_zone, NorS);

  // Let's create the maps of stratigraphy codes and K values, and Sc values  
  LSDLC.hardcode_fill_strati_to_K_map();
  LSDLC.hardcode_fill_strati_to_Sc_map();       

  // Now we'll get the lithology codes, convert them to K values and make a raster of K values for the model to use
  cout << "Getting the lithology codes" << endl;
  LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster);
  cout << "Converting the litho codes to K values" << endl;
  LSDRaster K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000001);
  cout << "Converting the litho codes to Sc values" << endl;
  LSDRaster Sc_values = LSDLC.index_to_Sc(lithocodes_index_raster, 0.45);
      
  LSDRaster Temp_raster = mod.return_as_raster();

  // Make a raster of uplift values 
  LSDRaster U_values;
  ifstream U_file_info_in;
  U_file_info_in.open(U_header.c_str());
  // check if the parameter file exists
  if( not U_file_info_in.fail() )
  {
    LSDRaster temp_U_values(U_fname,raster_ext);
    U_values = temp_U_values;
  }
  else
  {
    cout << "No variable U raster found. I will make one with a constant value of " << this_float_map["rudimentary_steady_forcing_uplift"] << endl;
    LSDRasterMaker URaster(Temp_raster);
    URaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_uplift"]);
    U_values = URaster.return_as_raster();
  }
  U_file_info_in.close();


  LSDRasterInfo RI(K_values);
  // Get the latitude and longitude
  cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
  LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

  int n_snaps;
  if (not this_bool_map["cycle_fluvial_snapping"])
  {
    n_snaps = 1;
  }
  else
  {
    n_snaps = this_int_map["snapping_cycles"];
  }

  int this_frame = 5000;
  mod.set_current_frame(this_frame);
  cout << "Let me do a few snaps. The number of snaps is: " << n_snaps << endl;


  LSDRaster lithocodes_raster;
  LSDRaster critical_slopes;
  LSDRaster new_muddpile_topo_raster;

  for (int i = 0; i < n_snaps; i++)
  {
    // Now the actual snapping!
    cout << endl << endl << endl << endl;
    cout << "=================================================" << endl;
    cout << "The U value at (500,500) is" << U_values.get_data_element(500,500) << endl;

    cout << "*****************************************" << endl;
    cout << "This is snap " << i+1 << " of " << n_snaps << endl;
    cout << "*****************************************" << endl;
    


    // if(i>0)
    // {
    cout << "Let me create a raster of the topography from the previous loop." << endl;
    new_muddpile_topo_raster = mod.return_as_raster();
    cout << "Now I'll check where in the lithocube we are and make a raster of lithocodes." << endl;
    lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(new_muddpile_topo_raster);
    cout << "I am now converting the lithocodes index raster to an LSDRaster" << endl;
    lithocodes_raster = LSDRaster(lithocodes_index_raster);
    cout << "Converting the litho codes to K values" << endl;
    K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000001);
    cout << "Converting the litho codes to Sc values" << endl;
    Sc_values = LSDLC.index_to_Sc(lithocodes_index_raster, 0.45);

    // }

    mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"]);

    // Snap to critical slopes now? 
    if (i == n_snaps-1)
    {
      // first, get the channels for burning
      LSDSpatialCSVReader chan_data = mod.get_channels_for_burning(this_int_map["threshold_contributing_pixels"]);

      cout << "Last snap! Getting critical slopes and smoothing" << endl;
      critical_slopes = mod.basic_valley_fill_critical_slope(Sc_values,this_int_map["threshold_contributing_pixels"]);



      // smooth the critical slope raster
      int kr = 0;
      float sigma = this_float_map["smoothing_window_radius"];
      for (int i = 0; i<this_int_map["smoothing_sweeps"]; i++)
      {
        cout << "smooth sweep " << i+1 << " of " << this_int_map["smoothing_sweeps"] << endl;
        critical_slopes = critical_slopes.GaussianFilter(sigma, kr);
      }

      // Now impose those channels
      LSDRasterMaker RM(critical_slopes);
      RM.impose_channels(chan_data);
      critical_slopes = RM.return_as_raster();

    }
    
    // Now print
    if (i % this_int_map["snap_print_rate"] == 0 || i == n_snaps-1)
    {
      cout << "This is snap " << i+1 << " of " << n_snaps << endl;

      if(this_bool_map["print_snapped_to_steady_frame"])
      {
        cout << endl << endl << "====================================" << endl;
        cout << "I am printing the snapped surface, increasing the frame count by 1." << endl;
        cout << "The frame number is " << this_frame << endl;
        mod.set_current_frame(this_frame);
        mod.print_rasters_and_csv( this_frame );

        if (i == n_snaps-1)
        {
          cout << "I'm writing the critical slope raster" << endl;
          string this_raster_name = OUT_DIR+OUT_ID+"_CritSlope";
          critical_slopes.write_raster(this_raster_name + to_string(this_frame), "bil");  

          cout << "Let me print the hillshade for you. " << endl;
          float hs_azimuth = 315;
          float hs_altitude = 45;
          float hs_z_factor = 1;
          LSDRaster hs_raster = critical_slopes.hillshade(hs_altitude,hs_azimuth,hs_z_factor);
          this_raster_name = OUT_DIR+OUT_ID+"_CritSlope_hs";
          hs_raster.write_raster(this_raster_name + to_string(this_frame), "bil");
        }

        if(this_bool_map["print_K_raster"])
        {
          cout << "Printing a K raster" << endl;
          K_values.write_raster(K_name + to_string(this_frame), "bil");
        }
        if(this_bool_map["print_Sc_raster"])
        {
          cout << "Printing a Sc raster" << endl;
          Sc_values.write_raster(Sc_name + to_string(this_frame), "bil");
        }
        if(this_bool_map["print_lithocode_raster"]) //&& i>0
        {
          cout << "Printing a lithocode raster" << endl;
          lithocodes_raster.write_raster(lithocodes_raster_name + to_string(this_frame), "bil");
        }

      }
    }
    this_frame++;

  }
}






























  if(this_bool_map["snap_to_critical_slopes"])
  {
    cout << "Let me snap the DEM to some critical slope values" << endl;
    cout << "Node that this happens after snapping to spatialy heterogeneous fluvial steady state. " << endl;
    cout << "It will superimpose the result on previous fluvial steady state calculations." << endl;

    // First load the two variable K and variable U maps
    string Sc_fname = DATA_DIR+this_string_map["Sc_raster_fname"];
    string Sc_header = DATA_DIR+this_string_map["Sc_raster_fname"]+".hdr";

    int n_snaps = this_int_map["snapping_hillslope_cycles"];

    if(this_bool_map["cycle_hillslope_snapping"])
    {
      int this_frame = 7000;
      mod.set_current_frame(this_frame);
      cout << "Let me do a few hillslope snaps. The number of snaps is: " << n_snaps << endl;
      LSDRaster critical_slopes;
      LSDRaster smoothed_raster;
      float pixel_weighting = 2;


      LSDRaster Temp_raster = mod.return_as_raster();
      LSDRasterInfo RI(Temp_raster);
      // Get the latitude and longitude
      cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
      LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );


      for (int i = 0; i < n_snaps; i++)
      {
        // Now the actual snapping!
        cout << endl << endl << endl << endl;
        cout << "=================================================" << endl;
        cout << "Snap cycle " << i << endl;
        
        if(this_bool_map["use_Sc_raster"])
        {
          ifstream Sc_file_info_in;
          Sc_file_info_in.open(Sc_header.c_str());
          // check if the parameter file exists
          if( not Sc_file_info_in.fail() )
          {
            LSDRaster temp_Sc_values(Sc_fname,raster_ext);
            critical_slopes = mod.basic_valley_fill_critical_slope(temp_Sc_values,this_int_map["threshold_contributing_pixels"]);
          }
          else
          {
            cout << "No variable Sc raster found. I will use a constant value of: " << this_float_map["rudimentary_critical_slope"] << endl;
            critical_slopes = mod.basic_valley_fill_critical_slope(this_float_map["rudimentary_critical_slope"],
                                                                    this_int_map["threshold_contributing_pixels"]);
          }
          Sc_file_info_in.close();
        }
        else
        {
          critical_slopes = mod.basic_valley_fill_critical_slope(this_float_map["rudimentary_critical_slope"],
                                                                this_int_map["threshold_contributing_pixels"]);
        }

        // Update the model and smooth
        mod.set_raster_data(critical_slopes);

        LSDRaster smoothed_raster = mod.basic_smooth(pixel_weighting);
        mod.set_raster_data(smoothed_raster);


        // Now print
        if (i % this_int_map["snap_hillslope_print_rate"] == 0)
        {
          cout << "This is hillslope snap " << i+1 << " of " << n_snaps << endl;

          if(this_bool_map["print_snapped_to_steady_frame"])
          {
            cout << endl << endl << "====================================" << endl;
            cout << "I am printing the snapped surface, increasing the frame count by 1." << endl;
            cout << "The frame number is " << this_frame << endl;
            mod.set_current_frame(this_frame);
            mod.print_rasters_and_csv( this_frame );
          }
        }
        this_frame++;

      }
    }
    else
    {
      cout << "I am snapping the hillslopes once. Find the data with CrtiSlope in the filename." << endl;
      LSDRaster critical_slopes;
      if(this_bool_map["use_Sc_raster"])
      {
        ifstream Sc_file_info_in;
        Sc_file_info_in.open(Sc_header.c_str());
        // check if the parameter file exists
        if( not Sc_file_info_in.fail() )
        {
          LSDRaster temp_Sc_values(Sc_fname,raster_ext);
          critical_slopes = mod.basic_valley_fill_critical_slope(temp_Sc_values,this_int_map["threshold_contributing_pixels"]);
        }
        else
        {
          cout << "No variable Sc raster found. I will use a constant value of: " << this_float_map["rudimentary_critical_slope"] << endl;
          critical_slopes = mod.basic_valley_fill_critical_slope(this_float_map["rudimentary_critical_slope"],
                                                                  this_int_map["threshold_contributing_pixels"]);
        }
        Sc_file_info_in.close();
      }
      else
      {
        critical_slopes = mod.basic_valley_fill_critical_slope(this_float_map["rudimentary_critical_slope"],
                                                              this_int_map["threshold_contributing_pixels"]);
      }

      cout << "I'm writing the critical slope raster" << endl;
      string this_raster_name = OUT_DIR+OUT_ID+"_CritSlope";
      critical_slopes.write_raster(this_raster_name,raster_ext);  

      cout << "Let me print the hillshade for you. " << endl;
      float hs_azimuth = 315;
      float hs_altitude = 45;
      float hs_z_factor = 1;
      LSDRaster hs_raster = critical_slopes.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

      this_raster_name = OUT_DIR+OUT_ID+"_CritSlope_hs";
      hs_raster.write_raster(this_raster_name,raster_ext); 
    }
    



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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
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
              int NRows = mod.get_NRows();
              int NCols = mod.get_NCols();
              LSDRasterMaker KRaster1(mod);
              KRaster1.resize_and_reset(NRows,NCols,this_float_map["DataResolution"],this_min_K);
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
          case 2:
            {
              cout << "Spatial variation in K. Simple linear contact horizontally across where S half of raster has a higher erodibility." << endl;
              int NRows = mod.get_NRows();
              int NCols = mod.get_NCols();
              cout << "NRows: " << NRows << " NCols: " << NCols << " DataRes: " << this_float_map["DataResolution"] << endl;
              LSDRasterMaker KRaster(mod);
              KRaster.resize_and_reset(NRows,NCols,this_float_map["DataResolution"],this_min_K);
              KRaster.increase_south_half_raster_values(this_int_map["increase_amt"]);
              this_K_raster = KRaster.return_as_raster();
            }
            break;
          default:
            {
              cout << "K variation. The options are 0 == random squares" << endl;
              cout << "  0 == random squares" << endl;
              cout << "  1 == sine waves (I lied, at the moment this doesn't work--SMM Sept 2017)." << endl;
              cout << "  2 == south half of raster has higher K (FJC Jul 2018)" << endl;
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
              cout << "Spatial variation in U. Case 1." << endl;
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
          case 2:
            {
                cout << "Spatial variation in U. Case 2." << endl;
                int NRows = mod.get_NRows();
                int NCols = mod.get_NCols();
                cout << "NRows: " << NRows << " NCols: " << NCols << " DataRes: " << this_float_map["DataResolution"] << endl;
                LSDRasterMaker URaster(mod);
                URaster.resize_and_reset(NRows,NCols,this_float_map["DataResolution"],this_float_map["min_U_for_spatial_var"]);
                URaster.increase_north_half_raster_values(this_int_map["increase_amt"]);
                this_U_raster = URaster.return_as_raster();
            }

            break;
          default:
            {
              cout << "Spatial variation in U. The options are 0 == random squares" << endl;
              cout << "  0 == random squares (I lied, at the moment this doesn't work--SMM Sept 2017)." << endl;
              cout << "  1 == sine waves" << endl;
              cout << "  2 == uplift top half of the raster at an order of magnitude faster than the bottom half (FJC June 2018)" << endl;
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
      int NRows = mod.get_NRows();
      int NCols = mod.get_NCols();
      cout << "NRows: " << NRows << " NCols: " << NCols << " DataRes: " << this_float_map["DataResolution"] << endl;
      LSDRasterMaker URaster(mod);
      URaster.resize_and_reset(NRows,NCols,this_float_map["DataResolution"],this_float_map["min_U_for_spatial_var"]);
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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
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
      cout << "Starting the variable uplift..." << endl;
      current_end_time = current_end_time+this_float_map["spatial_variation_time"];
      mod.set_endTime(current_end_time);

      mod.run_components_combined(this_U_raster, this_K_raster,use_adaptive_timestep);
      cout << "Now filling the raster" << endl;
      mod.raise_and_fill_raster();
    }

    if (this_bool_map["finish_with_steady_forcing"])   // FJC - added this because I want to run with a uniform uplift rate again afterwards but
    {                                                       // didn't want to mess up what we already had.
      cout << "Now running at steady state again" << endl;
      // now run at steady condition for a few extra cycles
      current_end_time = current_end_time+float(this_int_map["spatial_cycles"])*this_float_map["spatial_variation_time"];
      mod.set_uplift(0, this_float_map["min_U_for_spatial_var"]);
      mod.set_endTime(current_end_time);
      mod.run_components_combined();
    }
    else
    {
      // now run at steady condition for a few extra cycles
      current_end_time = current_end_time+float(this_int_map["spatial_cycles"])*this_float_map["spatial_variation_time"];
      mod.set_endTime(current_end_time);
      mod.run_components_combined(this_U_raster, this_K_raster,use_adaptive_timestep);
    }
  }

  //----------------------------------------------------------------------------------------------//
  // A simple fault uplift scenario where we just raise the north half of the raster by
  // a specified number of metres
  // FJC
  //----------------------------------------------------------------------------------------------//
  if (this_bool_map["simple_fault_uplift"])
  {
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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
    }

    float my_K=0.0001;
    // Calculate or set the K parameter depending on your choices about the simulation
    if(this_bool_map["calculate_K_from_relief"])
    {
      cout << "I am calculating a K value that will get a relief of " << this_float_map["fixed_relief"] << " metres" << endl;
      cout << " for an uplift rate of " << this_float_map["uplift_rate"]*1000 << " mm/yr" << endl;
      my_K = mod.fluvial_calculate_K_for_steady_state_relief(this_float_map["uplift_rate"],this_float_map["fixed_relief"]);
    }
    else
    {
      cout << "I am using a fixed K value of " << this_float_map["background_K"] << endl;
      my_K = this_float_map["background_K"];
    }

    // snap to steady state
    //mod.fluvial_snap_to_steady_state(this_float_map["uplift_rate"]);

    cout << "Starting the fault uplift. Run at steady forcing first" << endl;

    current_end_time = current_end_time+this_int_map["uplift_time"];

    mod.set_timeStep( this_float_map["dt"] );
    mod.set_K(my_K);

    cout << "The end time is " << current_end_time << endl;
    mod.set_endTime(current_end_time);

    mod.set_uplift(0, this_float_map["uplift_rate"]);
    mod.run_components_combined();

    cout << "Now adding the fault uplift. The north half of the raster will be uplifted by " << this_int_map["uplift_amt"] << " m." << endl;

    current_end_time = current_end_time+this_int_map["uplift_time"];
    mod.set_endTime(current_end_time);

    // normal fault across top two thirds of the raster.
    mod.normal_fault_part_of_raster(this_int_map["uplift_amt"], 2);
    mod.raise_and_fill_raster();

    mod.run_components_combined();
  }

  //----------------------------------------------------------------------------------------------//
  // A simple base level fall scenario
  // FJC
  //----------------------------------------------------------------------------------------------//
  if (this_bool_map["base_level_fall"])
  {
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
      if (this_bool_map["nonlinear"])
      {
        mod.set_nonlinear(true);
      }
    }

    float my_K=0.0001;
    // Calculate or set the K parameter depending on your choices about the simulation
    if(this_bool_map["calculate_K_from_relief"])
    {
      cout << "I am calculating a K value that will get a relief of " << this_float_map["fixed_relief"] << " metres" << endl;
      cout << " for an uplift rate of " << this_float_map["uplift_rate"]*1000 << " mm/yr" << endl;
      my_K = mod.fluvial_calculate_K_for_steady_state_relief(this_float_map["uplift_rate"],this_float_map["fixed_relief"]);
    }
    else
    {
      cout << "I am using a fixed K value of " << this_float_map["background_K"] << endl;
      my_K = this_float_map["background_K"];
    }

    // snap to steady state
    //mod.fluvial_snap_to_steady_state(this_float_map["uplift_rate"], this_float_map["fixed_relief"]);

    // cout << "Starting the base level fall scenario. Run at steady forcing first" << endl;
    //
    // current_end_time = current_end_time+this_int_map["uplift_time"];
    //
    mod.set_timeStep( this_float_map["dt"] );
    mod.set_K(my_K);
    //
    // cout << "The end time is " << current_end_time << endl;
    // mod.set_endTime(current_end_time);
    //
    // mod.set_uplift(0, this_float_map["uplift_rate"]);
    // mod.run_components_combined();

    cout << "Now dropping the base level by " << this_int_map["uplift_amt"] << " m." << endl;

    // reduce the print interval
    mod.set_print_interval( this_int_map["transient_print_interval"] );

    current_end_time = current_end_time+this_int_map["uplift_time"];
    mod.set_endTime(current_end_time);
    //mod.raise_and_fill_raster(this_float_map["min_slope_for_fill"]);

    mod.base_level_fall(this_int_map["uplift_amt"]);
    mod.raise_and_fill_raster(this_float_map["min_slope_for_fill"]);

    mod.run_components_combined();
  }

//--------------------------------------------------------------------------------------//
// Run the tilting scenario
// FJC 11/03/19
//--------------------------------------------------------------------------------------//
if( this_bool_map["run_progressive_tilt"] )
{
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
    if (this_bool_map["nonlinear"])
    {
      mod.set_nonlinear(true);
    }
  }

  // set boundary conditions based on the tilt boundary.
  vector <string> bc(4, "n");          // Initialise boundaries to No flow
  if (this_string_map["tilt_boundary"] == "N")
  {
    bc[0] = "b";
    bc[1] = "p";
    bc[2] = "n";
    bc[3] = "p";
  }
  else if (this_string_map["tilt_boundary"] == "E")
  {
    bc[0] = "p";
    bc[1] = "b";
    bc[2] = "p";
    bc[3] = "n";
  }
  else if (this_string_map["tilt_boundary"] == "S")
  {
    bc[0] = "n";
    bc[1] = "p";
    bc[2] = "b";
    bc[3] = "p";
  }
  else if (this_string_map["tilt_boundary"] == "W")
  {
    bc[0] = "p";
    bc[1] = "n";
    bc[2] = "p";
    bc[3] = "b";
  }
  mod.set_boundary_conditions( bc );          // Set these as default boundary conditions

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

  // run progressive tilting for a defined period of time (set uplift field to tilted block)
  cout << "Running progressive tilting of " << this_float_map["tilt_angle"] << " deg, for " << this_int_map["tilt_time"] << " years" << endl;

  // make a U raster based on a tilted block
  LSDRasterMaker URaster(this_int_map["NRows"],this_int_map["NCols"]);
  URaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],this_float_map["snapped_to_steady_uplift"]);
  URaster.tilted_block(this_float_map["tilt_angle"], this_string_map["tilt_boundary"]);
  LSDRaster this_U_raster = URaster.return_as_raster();
  // write the raster
  string U_fname = OUT_DIR+OUT_ID+"_URaster";
  string bil_name = "bil";
  this_U_raster.write_raster(U_fname,bil_name);

  // make a spatially invariant K raster.
  LSDRasterMaker KRaster(this_int_map["NRows"],this_int_map["NCols"]);
  KRaster.resize_and_reset(this_int_map["NRows"],this_int_map["NCols"],this_float_map["DataResolution"],new_K);
  LSDRaster this_K_raster = KRaster.return_as_raster();

  // reduce the print interval
  mod.set_print_interval( this_int_map["tilt_print_interval"] );

  current_end_time = current_end_time+this_int_map["tilt_time"];
  mod.set_endTime(current_end_time);

  mod.run_components_combined(this_U_raster, this_K_raster, this_bool_map["use_adaptive_timestep"]);

  // now run with a spatially invariant uplift block until steady state
  current_end_time = current_end_time+this_int_map["uplift_time"];
  mod.set_endTime(current_end_time);
  mod.set_uplift(0, this_float_map["snapped_to_steady_uplift"]);
  mod.run_components_combined();
}

if( this_bool_map["run_instantaneous_tilt"] )
{
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
    if (this_bool_map["nonlinear"])
    {
      mod.set_nonlinear(true);
    }
  }

  // set boundary conditions based on the tilt boundary.
  vector <string> bc(4, "n");          // Initialise boundaries to No flow
  if (this_string_map["tilt_boundary"] == "N")
  {
    bc[0] = "b";
    bc[1] = "p";
    bc[2] = "n";
    bc[3] = "p";
  }
  else if (this_string_map["tilt_boundary"] == "E")
  {
    bc[0] = "p";
    bc[1] = "b";
    bc[2] = "p";
    bc[3] = "n";
  }
  else if (this_string_map["tilt_boundary"] == "S")
  {
    bc[0] = "n";
    bc[1] = "p";
    bc[2] = "b";
    bc[3] = "p";
  }
  else if (this_string_map["tilt_boundary"] == "W")
  {
    bc[0] = "p";
    bc[1] = "n";
    bc[2] = "p";
    bc[3] = "b";
  }
  mod.set_boundary_conditions( bc );          // Set these as default boundary conditions

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

  cout << "Now I'm tilting your landscape by " << this_float_map["tilt_angle"] << " deg" << endl;

  // reduce the print interval
  mod.set_print_interval( this_int_map["tilt_print_interval"] );

  current_end_time = current_end_time+this_int_map["uplift_time"];
  mod.set_endTime(current_end_time);
  //mod.raise_and_fill_raster(this_float_map["min_slope_for_fill"]);

  mod.instantaneous_tilt(this_float_map["tilt_angle"], this_string_map["tilt_boundary"]);
  mod.run_components_combined();
  }
}
