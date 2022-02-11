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
#include "../LSDChiTools.hpp"
#include "../LSDIndexRaster.hpp"
using namespace std;


int main (int nNumberofArgs,char *argv[])
{

  string version_number = "0.08d";
  string citation = "https://www.doi.org/10.5281/zenodo.997388";


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
  cout << "|| This is MuddPILE version                                   ||" << endl;
  cout << "|| " << version_number  << endl;
  cout << "|| If the version number has a d at the end it is a           ||" << endl;
  cout << "||  development version.                                      ||" << endl;
  cout << "================================================================" << endl;


  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);
  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // Check if we are doing the version or the citation
  if(f_name == "lsdtt_citation.txt")
  {

    cout << endl << endl << endl << "==============================================" << endl;
    cout << "To cite this code, please use this citation: " << endl;
    cout << citation << endl;
    cout << "Copy this url to find the full citation." << endl;
    cout << "also see above for more detailed citation information." << endl;
    cout << "=========================================================" << endl;

    ofstream ofs;
    ofs.open("./muddpiledriver-citation.txt");
    ofs << citation << endl;
    ofs.close();

    exit(0);
  }

  if(f_name == "lsdtt_version.txt")
  {
    ofstream ofs;
    ofs.open("./muddpiledriver-version.txt");
    ofs << version_number << endl;
    ofs.close();

    exit(0);
  }

  // load parameter parser object
  LSDParameterParser LSDPP(path_name,f_name);

  // for the chi tools we need georeferencing so make sure we are using bil format
  LSDPP.force_bil_extension();

  // maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

  // this will contain the help file
  map< string, vector<string> > help_map;


  //==================================================================================
  //
  // .#####....####...#####....####...##...##..######..######..######..#####....####..
  // .##..##..##..##..##..##..##..##..###.###..##........##....##......##..##..##.....
  // .#####...######..#####...######..##.#.##..####......##....####....#####....####..
  // .##......##..##..##..##..##..##..##...##..##........##....##......##..##......##.
  // .##......##..##..##..##..##..##..##...##..######....##....######..##..##...####..
  //
  //=================================================================================
  // Parameters for initiating the model domain
  int_default_map["NRows"] = 200;
  help_map["NRows"] = {  "int","200","Number of rows in a model where you have not loaded a DEM.","This is overwritten if you load an initial DEM."};

  int_default_map["NCols"] = 400;
  help_map["NCols"] = {  "int","400","Number of cols in a model where you have not loaded a DEM.","This is overwritten if you load an initial DEM."};

  float_default_map["DataResolution"] = 30;
  help_map["DataResolution"] = {  "float","30","Pixel resolution of a model where you have not loaded a DEM.","This is overwritten if you load an initial DEM."};

  // Parameters that are used if you load a raster
  bool_default_map["read_initial_raster"] = false;
  help_map["read_initial_raster"] = {  "bool","false","Reads a raster from a DEM using the read path and read fname.","Make sure your raster in in ENVI bi format and in UTM coordinate system EPSG:326XX or EPSG:327XX."};
 
  float_default_map["minimum_elevation"] = 0.0;
  help_map["minimum_elevation"] = { "float","0.0","All elevation values below this become nodata if remove_seas is true.","Ususally 0."};

  float_default_map["maximum_elevation"] = 30000;
  help_map["maximum_elevation"] = {  "float","0.0","All elevation values above this become nodata if remove_seas is true.","Pick a big number."};

  float_default_map["min_slope_for_fill"] = 0.0001;
  help_map["min_slope_for_fill"] = {  "float","0.0001","Minimum slope between pixels for the filling algorithm.","Best not to change the default."};

  bool_default_map["remove_seas"] = true; // elevations above minimum and maximum will be changed to nodata
  help_map["remove_seas"] = {  "bool","true","Slightly misleading name; it replaces both high and low DEM values with nodata.","This gets rid of low lying areas but also is handy when the nodata is not translated from the raw DEM and it is full of funny large numbers."};
 
  bool_default_map["print_raster_without_seas"] = false;
  help_map["print_raster_without_seas"] = {  "bool","false","Overwrites the raster without seas.","DANGER this will replace your existing raster with the high and low points replaced by nodata. See the remove_seas flag"};

  bool_default_map["convert_initial_raster_to_diamond_square_fractal"] = false;
  help_map["convert_initial_raster_to_diamond_square_fractal"] = {  "bool","false","This takes the initial raster (with its nodata footprint) and creates a random diamond square raster on to original footprint.","The spinup of the diamond square surface will follow the diamond square spinup parameters elsewhere in this program."};

  bool_default_map["test_single_basin_outlet"] = false;
  help_map["test_single_basin_outlet"] = {  "bool","false","This tests the buffering routines around the outside of the DEM.","Most useful for the diamond square conversion."};

  bool_default_map["raise_raster_over_base_level"] = true;
  help_map["raise_raster_over_base_level"] = {  "bool","false","This looks at the baselevel file and if the lowest point is above the lowest point in the model raster then the raster is raised above the baselevel.","If pixels are above the base level channel then the raster will not incise."};



  // Some logic for filling rasters that are being loaded
  bool_default_map["carve_before_fill"] = false; // Implements a carving algorithm
  help_map["carve_before_fill"] = {  "bool","false","This implements a breaching algorithm before filling.","Good for landscapes with DEM obstructions (like roads) across the channels."};

  bool_default_map["raster_is_filled"] = false;
  help_map["raster_is_filled"] = {  "bool","false","This reads a pre-existing fill raster to save time.","You need to have printed the fill raster if you set this to true."};

  bool_default_map["print_fill_raster"] = false;
  help_map["print_fill_raster"] = {  "bool","false","Prints the fill raster.","Filename includes _FILL"};

 // Read a channel file
  bool_default_map["extract_single_channel"] = false;
  help_map["extract_single_channel"] = {  "bool","false","Extracts a single channel from the DEM. The channel is a flow path that starts from a point designated by channel_source_fname.","Doing this in lsdtt-basic-metrics provides more options."};

  string_default_map["channel_source_fname"] = "channel_source_fname";
  help_map["channel_source_fname"] = {  "string","channel_source_fname","A csv with latitude and longitude as headers that is the source point of a flow path that will be printed if extract_single_channel is true.","Only reads the first point so far."};

  string_default_map["single_channel_print_fname"] = "single_channel";
  help_map["single_channel_print_fname"] = {  "bool","single_channel","The file prefix of the single channel csv. This has data from a single flow path.","Columns include latitude longitude flow distance, area, elevation and a few other items"};

  // paramters for controlling model output
  int_default_map["print_interval"] = 10;
  help_map["print_interval"] = {  "int","10","The number of timesteps required for printing.","This is reactive to the adaptive timestep so will give a fixed printing interval."};
  
  bool_default_map["write_hillshade"] = true;
  help_map["write_hillshade"] = {  "bool","false","Write the hillshade raster.","You need this for a lot of our plotting routines. Filename includes _HS"};


  // control of m and n, and paramters for chi
  float_default_map["A_0"] = 1;
  help_map["A_0"] = {  "float","1.0","The A_0 parameter for chi computation.","Usually set to 1 so that the slope in chi-elevation space is the same as k_sn"};
 
  float_default_map["m"] = 0.5;
  help_map["m"] = {  "float","0.5","The area exponent in the stream power law.","Defaults lead to m/n = 0.5"};

  float_default_map["n"] = 1;
  help_map["n"] = {  "float","1","The slope exponent in the stream power law.","Model slows down a lot if n does not equal 1 but there is much evidence that n is greater than 1 for example Harel et al 2016 Geomorphology."};

  int_default_map["uplift_mode"] = 0;
  help_map["uplift_mode"] = {  "int","0","Model can run with specific uplift fields. This is a bit old and we now use uplift rasters.","default == block; uplift 1 == tilt block; 2 == gaussian; 3 == quadratic; 4 == periodic."};

  float_default_map["dt"] = 250;
  help_map["dt"] = {  "float","250","The starting timestep in years.","This can change if you use adaptive timesteps."};

  float_default_map["D"] = 0.002;
  help_map["D"] = {  "float","0.002","Default sediment transport coefficient.","In m^2 per year"};

  float_default_map["S_c"] = 1.0;
  help_map["S_c"] = {  "float","1.0","Default critical slope parameter for hillslopes.","Dimensionless"}; 

  float_default_map["background_K"] = 0.0005;
  help_map["background_K"] = {  "float","0.0005","Default channel erodibility.","Dimensions depend on m and n"}; 

  // Parameters for the initial surface
  bool_default_map["use_diamond_square_initial"] = true;
  help_map["use_diamond_square_initial"] = {  "bool","true","Uses a diamond square algorithm to produce a pseudo fractal surface as the initial surface in a rectangular model run.","You get better branching behaviour using this algorithm in the channel network than if you use other initial surfaces. You can now use diamond_square_spinup and this is preferred"}; 

  float_default_map["diamond_square_relief"] = 16;
  help_map["diamond_square_relief"] = {  "float","16","The total releif of the original diamond square surface after initiation.","Usually you snap to stead after this to control the final relief but this controls how oriented the initial drainage patterns are."}; 

  int_default_map["diamond_square_feature_order"] = 8;
  help_map["diamond_square_feature_order"] = {  "int","8","Diamond square has fractal features that are of order 2^n and this sets the maximum value of n.","You will want 2^n to be less than the number of rows or columns (whichever is less) in the raster. Sets the repeating scale of the big features."}; 

  int_default_map["ds_minimum_feature_order"] = 4;
  help_map["ds_minimum_feature_order"] = {  "int","4","Sets the minimum scale of the diamon square features. This is 2^n where n is the feature order.","I wouldn't change this."}; 

  bool_default_map["taper_edges"] = true;
  help_map["taper_edges"] = {  "bool","true","This reduces the relief along the edges so you don't get weird negative elevations below base level along the edge.","Setting to false might lead to some weird behaviour so don't change."}; 

  int_default_map["taper_rows"] = 10;
  help_map["taper_rows"] = {  "int","10","Number of rows over which the model tapers the initial surface to 0 at the edge of the DEM","Larger numbers mean smoother edges."}; 


  bool_default_map["superimpose_parabola"] = true;
  help_map["superimpose_parabola"] = {  "bool","true","After initial construction of the diamond square surface this adds a parabolic surface to get channel to drain away from the centre","If this is false you will get a lot of lakes in the middle of the initial DEM."}; 

  float_default_map["parabola_relief"] = 6;
  help_map["parabola_relief"] = {  "float","6","Relief of the parabolic surface added to the DEM","Too small and you get a lot of lakes. Too big and the channels start to straighten our toward the edge."}; 

  bool_default_map["roughen_surface"] = true;
  help_map["roughen_surface"] = {  "bool","true","Adds some noise after the initial DEM", "This tries to add some roughness to fill steps"}; 

  bool_default_map["fill_then_roughen_surface"] = true;
  help_map["fill_then_roughen_surface"] = {  "bool","true","The initial diamond square surface will produce some low topography that will get filled. Adter filling you will get very straight channels. So this roughens the surface before filling to get a more random network in the filled lakes", "Used to deal with filled depressions after inital diamond square step."}; 

  float_default_map["roughness_relief"] = 0.25;
  help_map["roughness_relief"] = {  "float","0.25","Amplitude of diamond square roughening elements. Uses a unifrom distribution", "Used to deal with filled depressions after inital diamond square step."}; 

  bool_default_map["diffuse_initial_surface"] = false;
  help_map["diffuse_initial_surface"] = {  "bool","false","Smooths initial surface", "Not sure why I had this--SMM"}; 

  int_default_map["diffuse_steps"] = 5;
  help_map["diffuse_steps"] = {  "int","5","Number of smoothing cycles", "More means smoother"}; 

  bool_default_map["print_initial_surface"] = true;
  help_map["print_initial_surface"] = {  "bool","true","Does what is says on the tin", "Self explanitory"}; 


  // Parameters for spinning up the simulation
  bool_default_map["spinup"] = false;
  help_map["spinup"] = {  "bool","false","For spinning up a rectangular model.","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["spinup_K"] = 0.001;
  help_map["spinup_K"] = {  "float","0.001","For spinning up a rectangular model, the K value.","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["spinup_U"] = 0.001;
  help_map["spinup_U"] = {  "float","0.001","For spinning up a rectangular model, the U value.","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["spinup_dt"] = 250;
  help_map["spinup_dt"] = {  "float","250","For spinning up a rectangular model, the dt value.","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["spinup_time"] = 20000;
  help_map["spinup_time"] = {  "float","20000","For spinning up the time to spend spinning up the model.","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  bool_default_map["staged_spinup"] = true;
  help_map["spinup"] = {  "bool","true","For spinning up a rectangular model. I forgot what this does but I don't use it anymore-SMM","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 


  // This spinup routine tries to combine snapping and hillslope diffusion
  bool_default_map["cyclic_spinup"] = false;
  help_map["cyclic_spinup"] = {  "bool","false","For spinning up a rectangular model. Turns hillslopes on and off sequentially","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  int_default_map["spinup_cycles"] = 5;
  help_map["spinup_cycles"] = {  "int","5","For spinning up a rectangular model. Turns hillslopes on and off sequentially. This is the number of cycles","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  bool_default_map["cycle_K"] = true;
  help_map["cycle_K"] = {  "bool","true","For spinning up a rectangular model. This varies K in the cycles","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  bool_default_map["cycle_U"] = false;
  help_map["cycle_U"] = {  "bool","false","For spinning up a rectangular model. This varies U in the cycles","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["cycle_K_factor"] = 2;
  help_map["cycle_K_factor"] = {  "float","2","For spinning up a rectangular model. This factor by which K varies in the cycles","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 

  float_default_map["cycle_U_factor"] = 2;
  help_map["cycle_U_factor"] = {  "float","2","For spinning up a rectangular model. This factor by which U varies in the cycles","More useful to use diamond_square_spinup and snapping so not advised to use this"}; 



  // The diamond square spinup routine
  // This generates the nicest initial surfaces.
  bool_default_map["diamond_square_spinup"] = false;
  help_map["diamond_square_spinup"] = {  "bool","false","For spinup of a rectangular model run. Uses the diamond square algorithm as an initial surface.","This is done so the channel network is dendritic"};

  // control of snapping to steady state
  bool_default_map["snap_to_steady"] = false;
  help_map["snap_to_steady"] = {  "bool","false","This calculates the relief in the fluvial nework by solving stream power for a given uplift and K.","Used to start model runs with a steady landscape"};

  float_default_map["snapped_to_steady_uplift"] = 0.0001;
  help_map["snapped_to_steady_uplift"] = {  "float","0.0001","The uplift rate used for the snap_to_steady routine.","In  m/yr"};

  float_default_map["snapped_to_steady_relief"] = 400;
  help_map["snapped_to_steady_relief"] = {  "float","400","For the snap_to_steady option this back calculates K to get this relief in the landscape.","In  m"};

  bool_default_map["print_snapped_to_steady_frame"] = false;
  help_map["print_snapped_to_steady_frame"] = {  "bool","false","Prints the snapped DEM.","works with snap_to_steady"};

  // More complex snapping
  bool_default_map["snap_to_steady_variable_variable_K_variable_U"] = false;
  help_map["snap_to_steady_variable_variable_K_variable_U"] = {  "bool","false","A snapping routine with spatially variable U and K.","Assign K and U with rasters using variable_K_name and variable_U_name"};

  string_default_map["variable_K_name"] = "variable_K";  
  help_map["variable_K_name"] = {  "string","variable_K","For snap_to_steady_variable_variable_K_variable_U this is the name of the K raster.","There are some flags to do this with a constant K as well. Use make_constant_K_raster"};
 
  string_default_map["variable_U_name"] = "variable_U";
  help_map["variable_U_name"] = {  "string","variable_U","For snap_to_steady_variable_variable_K_variable_U this is the name of the U raster.","There are some flags to do this with a constant U as well. Use make_constant_U_raster"};
 
  // Even more complex snapping
  bool_default_map["snap_to_steady_variable_variable_K_variable_U_use_LithoCube"] = false;
  help_map["snap_to_steady_variable_variable_K_variable_U_use_LithoCube"] = {  "bool","false","A snapping routine with spatially variable U and K and uses a lithocube.","U raster is assigned and K raster is determined with lithocube parameters"};

  string_default_map["fixed_channel_csv_name"] = "NULL";  // This contains the fixed channel data
  help_map["fixed_channel_csv_name"] = {  "string","NULL","This is the name of a csv file with a fixed channel.","Channel needs column headers latitude longitude elevation and has csv in the name."};

  bool_default_map["only_test_single_channel"] = false;
  help_map["only_test_single_channel"] = {  "bool","false","This checks the format of the fixed channel and exits.","For bug checking."};

  bool_default_map["only_test_impose_single_channel"] = false;
  help_map["only_test_impose_single_channel"] = {  "bool","false","Imposes single channel prints DEM and exits.","For bug checking."};

  bool_default_map["buffer_single_channel"] = true;
  help_map["buffer_single_channel"] = {  "bool","true","Adds data to nodata pixels next to the single channel so that after flow routing the channel does not drain off the edge of the DEM.","Useful for when the single channel runs along the edge of the DEM."};

  bool_default_map["fixed_channel_dig_and_slope"] = false;
  help_map["fixed_channel_dig_and_slope"] = {  "bool","false","For the single channel this enforces the slope and digs the channel in a bit to ensure flow is routed through it","The digging is used to enforce flow routing"};

  float_default_map["single_channel_drop"] = 0;
  help_map["single_channel_drop"] = {  "float","0","For the single channel this drops the elevation","Use if you want a big drop in base level"};

  string_default_map["single_channel_fd_string"] = "flowdistance(m)";
  help_map["single_channel_fd_string"] = {  "string","flow distance(m)","The name of the flow distance column in the single channel csv.","Case sensitive"};

  string_default_map["single_channel_elev_string"] = "elevation(m)";
  help_map["single_channel_elev_string"] = {  "string","elevation(m)","The name of the elevation column in the single channel csv.","Case sensitive"};

  string_default_map["test_single_channel_name"] = "test_single_channel";
  help_map["test_single_channel_name"] = {  "string","test_single_channel","This is the filename of the single channel that prints after the slope and dig are imposed","File is used for bug checking."};

  bool_default_map["lower_raster_before_channel_buffer"] = false;
  help_map["lower_raster_before_channel_buffer"] = {  "bool","false","Lowers the raster before the single channel buffering happens","Used in case raster is dropped."};
 
  float_default_map["lower_raster_before_buffer_amount"] = 0;
  help_map["lower_raster_before_buffer_amount"] = {  "float","0","Lowers the raster before the single channel by this amount","Used in case raster is dropped."};
   
  bool_default_map["print_updated_fixed_channel"] = false;
  help_map["print_updated_fixed_channel"] = {  "bool","false","Prints the single channel after imposition of the fixed slope and diggin.","For bug checking."};


  // Yet more complex snapping that attempts to migrate divides by cycling snapping
  bool_default_map["cycle_fluvial_snapping"] = false;
  help_map["cycle_fluvial_snapping"] = {  "bool","false","This cycles the snapping so the divides can migrate in the event of spatially heterogeneuous lithology.","Allow divides to move."};

  int_default_map["snapping_cycles"] = 30;
  help_map["cycle_fluvial_snapping"] = {  "int","30","Number of fluvial snapping cycles.","Turn on with cycle_fluvial_snapping."};

  int_default_map["snap_print_rate"] = 10;     // how frequently the snapping is printed in terms of snap cycles
  help_map["snap_print_rate"] = {  "int","10","This prints the DEM evey n cycles set by this parameters.","Turn on with cycle_fluvial_snapping."};


  // Yet more complex snapping that attempts to migrate divides by cycling snapping of the critical slope component
  bool_default_map["cycle_hillslope_snapping"] = false;
  help_map["cycle_hillslope_snapping"] = {  "bool","false","This includes the hillslopes in the cyclic snapping.","Turn on with cycle_fluvial_snapping."};

  int_default_map["snapping_hillslope_cycles"] = 1;
  help_map["snapping_hillslope_cycles"] = {  "int","1","Not used in the lithosnap but in only hillslope snapping this is the number of snaps.","More for bug checking has no affect on lithosnap."};

  int_default_map["snap_hillslope_print_rate"] = 10;     // how frequently the snapping is printed in terms of snap cycles  
  help_map["snap_hillslope_print_rate"] = {  "int","10","Not used in the lithosnap but in only hillslope snapping this is printing for hillslope snaps.","More for bug checking has no affect on lithosnap."};


  // Snapping if you use the initial raster as a lid so that
  // the snapped landscape cannot go above the current surface
  bool_default_map["use_initial_topography_as_max_elevation"] = false;
  help_map["use_initial_topography_as_max_elevation"] = {  "bool","false","For snapping this caps the elevation at the elevation of the input DEM.","Used to make sure you don't snap above present elevations."};

  // This is for printing of the channel network
  bool_default_map["print_initial_channel_network"] = false;
  help_map["use_initial_topography_as_max_elevation"] = {  "bool","false","Prints the initial channel network after the base level has been imposed.","Used to bug check the base level snapping."};

  bool_default_map["print_final_channel_network"] = false;
  help_map["print_final_channel_network"] = {  "bool","false","Prints the final channel network after the base level has been imposed.","Used to bug check the base level snapping."};


  bool_default_map["impose_single_channel_printing"] = false;
  help_map["impose_single_channel_printing"] = {  "bool","false","Prints the imposed channel to csv. This is after it has been snapped to pixels and slope has been enforced.","Used to bug check the base level snapping."};

  bool_default_map["nodata_single_channel_printing"] = false;
  help_map["nodata_single_channel_printing"] = {  "bool","false","Prints the imposed channel to csv. Includes nodata to identify problem pixels.","Used to bug check the base level snapping."};

  bool_default_map["convert_csv_to_geojson"] = false;
  help_map["convert_csv_to_geojson"] = {  "bool","false","Converts csv files to geojson files","Makes csv output easier to read with a GIS. Warning: these files are much bigger than csv files."};

  // This is for the critical slope
  float_default_map["rudimentary_critical_slope"] = 0.8;
  help_map["rudimentary_critical_slope"] = {  "float","0.8","If you don't want to use the lithocube for critical slopes this just sets S_c the same everywhere to this value.","Normally you would use the lithocube but this is for testing routines without needing the lithocube which takes a long time to load."};

  int_default_map["threshold_contributing_pixels"] = 1000;
  help_map["threshold_contributing_pixels"] = {  "int","1000","Threshold number of accumulated pixels before you get a hillslope.","Pixels not area"};

  bool_default_map["snap_to_critical_slopes"] = false;
  help_map["snap_to_critical_slopes"] = {  "bool","false","A bug testing routine that snaps to critical slopes but does not fluvial snapping.","For bug testing"};

  string_default_map["Sc_raster_fname"] = "S_c";
  help_map["Sc_raster_fname"] = {  "string","S_c","If you load the critical slope raster from file this is its prefix (without bil extension).","You can combine this with creating an S_c raster"};

  bool_default_map["use_Sc_raster"] = false;
  help_map["use_Sc_raster"] = {  "bool","false","Loads an S_c raster for critical slopes rather than computing from the rudimentary value or from the lithocube","You can combine this with creating an S_c raster"};
 
  bool_default_map["snap_to_steady_critical_slopes_use_LithoCube"] = false;
  help_map["snap_to_steady_critical_slopes_use_LithoCube"] = {  "bool","false","A bug testing routine that snaps to critical slopes but does not fluvial snapping. THis one tests the lithocube","For bug testing"};

  bool_default_map["print_Sc_raster"] = false; 
  help_map["use_Sc_raster"] = {  "bool","false","Prints S_c raster used in a run","For checking your S_c values"};
 
  string_default_map["temp_chan_name_for_crit_slopes"] = "temp_channels";
  help_map["temp_chan_name_for_crit_slopes"] = {  "string","temp_channels","Due to some nuances of the code you must write the channel to file and then read the file for the hillslope snapping. If you are batch processing you will need to have different names for the channel files","Batch processing scripts should make sure this name is different for each model implementation."};
 
  // Smoothing parameters 
  float_default_map["smoothing_window_radius"] = 50;
  help_map["smoothing_window_radius"] = {  "float","50","For snapping this smooths the raster at the end of the snapping cycles over a window of this radius.","In metres"};

  int_default_map["smoothing_sweeps"] = 1;
  help_map["smoothing_sweeps"] = {  "int","1","The number of smoothing iterations if you smooth after snapping.","More diffused landscape if this number is higher"};

  // forcing of dissection
  bool_default_map["force_dissect"] = true;
  help_map["force_dissect"] = {  "bool","true","For spinup options and diamond square this forces channel dissection.","Used to ensure channel connectivity"};

  int_default_map["force_dissect_steps"] = 10000;
  help_map["force_dissect_steps"] = {  "int","10000","For spinup options and diamond square this forces the number channel dissection iterations. The higher the number the better your chance of getting a fully dissected landscape.","Used to ensure channel connectivity"};

  // Some parameters for very rudimentary steady forcing
  bool_default_map["make_constant_K_raster"] = false;
  help_map["make_constant_K_raster"] = {  "bool","false","Makes a K raster that can be read by model runs.","The value of K is set by rudimentary_steady_forcing_K"};

  bool_default_map["make_constant_U_raster"] = false;
  help_map["make_constant_U_raster"] = {  "bool","false","Makes a U raster that can be read by model runs.","The value of K is set by rudimentary_steady_forcing_uplift"};

  bool_default_map["make_constant_Sc_raster"] = false;
  help_map["make_constant_Sc_raster"] = {  "bool","false","Makes a Sc raster that can be read by model runs.","The value of Sc is set by rudimentary_critical_slope"};

  bool_default_map["rudimentary_steady_forcing"] = false;
  help_map["rudimentary_steady_forcing"] = {  "bool","false","This just runs a single forcing for a fixed time.","See other default parameters with rudimentary in the name."};

  float_default_map["rudimentary_steady_forcing_time"] = 100000;
  help_map["rudimentary_steady_forcing_time"] = {  "float","100000","The time of a rudimentary forced run.","See other default parameters with rudimentary in the name."};

  float_default_map["rudimentary_steady_forcing_uplift"] = 0.0005;
  help_map["rudimentary_steady_forcing_uplift"] = {  "float","0.0005","Uplift rate in the rudimentary forced runs and also used as a fixed uplift field.","Units are m/yr."};

  float_default_map["rudimentary_steady_forcing_K"] = 0.000005;
  help_map["rudimentary_steady_forcing_K"] = {  "float","0.000005","The K value used in the constant K raster.","Units depend on m and n."};

  // Parameters for hillslopes
  bool_default_map["hillslopes_on"] = false;
  help_map["hillslopes_on"] = {  "bool","false","Turn the hillslope module on.","This is computationally intensive and doesn't work on irregular boundaries."};

  bool_default_map["nonlinear"] = false;
  help_map["nonlinear"] = {  "bool","false","Use the nonlinear hillslope module.","This is the roering et al model."};

  bool_default_map["use_hillslope_hybrid"] = false;
  help_map["use_hillslope_hybrid"] = {  "bool","false","If true the model will use a hybrid hillslope model.","This is linear diffusion combined with a critical slope."};



  // BELOW ARE SETTINGS FOR DIFFERENT RECTANGULAR MODEL SCENARIOS
  // Some parameters for cyclic forcing
  // these also inherit parameters from the cyclic spinup
  bool_default_map["run_cyclic_forcing"] = false;
  help_map["run_cyclic_forcing"] = {  "bool","false","For rectangular model runs. Cycles the forcing (either U or K) through time.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};

  float_default_map["cyclic_forcing_time"] = 10000;
  help_map["cyclic_forcing_time"] = {  "float","1000","The temporal wavelength of cyclic forcing.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};

  float_default_map["baseline_K_for_cyclic"] = 0.0001;
  help_map["baseline_K_for_cyclic"] = {  "float","0.0001","The lowest K value in the cyclic K forcing. Maximum determined by cycle_K_factor.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};

  float_default_map["baseline_U_for_cyclic"] = 0.0005;
  help_map["baseline_U_for_cyclic"] = {  "float","0.0005","The lowest U value in the cyclic U forcing. Maximum determined by cycle_UY_factor.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};

  int_default_map["cyclic_cycles"] = 3;
  help_map["cyclic_cycles"] = {  "int","3","Number of cycles in the cyclic forcing.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};

  float_default_map["cyclic_dt"] = 250;
  help_map["cyclic_dt"] = {  "float","250","Timestep during cyclic forcing.","Forcings follow a sine wave. Examples in Mudd 2017 ESPL."};


  // some parameters for setting K to a fixed uplift and relief
  bool_default_map["set_fixed_relief"] = false;
  help_map["set_fixed_relief"] = {  "bool","false","Back calculates K to force a value of relief in the initial model.","Used to make sure different model domains have similar relief."};


  float_default_map["fixed_relief"] = 500;
  help_map["fixed_relief"] = {  "float","500","Set the relief of the steady state run.","Once relief if set this back calculates K."};



  // some parameters for having a random uplift through time
  bool_default_map["run_random_forcing"] = false;
  help_map["run_random_forcing"] = {  "bool","false","Runs a transient scenario with random periods of uplift","For running models with random uplift histories"};

  float_default_map["maximum_time_for_random_cycle"] = 20000;
  help_map["maximum_time_for_random_cycle"] = {  "float","20000","For random uplift this is the maximum time of a step","For running models with random uplift histories"};

  float_default_map["minimum_time_for_random_cycle"] = 5000;
  help_map["minimum_time_for_random_cycle"] = {  "float","5000","For random uplift this is the minimum time of a step","For running models with random uplift histories"};

  float_default_map["maximum_U_for_random_cycle"] = 0.001;
  help_map["maxcimum_U_for_random_cycle"] = {  "float","0.001","For transient run with random uplift this is the maximum uplift.","For running models with random uplift histories."};

  float_default_map["minimum_U_for_random_cycle"] = 0.0001;
  help_map["minimum_U_for_random_cycle"] = {  "float","0.0001","Set the relief of the steady state run.","For running models with random uplift histories"};

  float_default_map["random_dt"] = 10;
  help_map["random_dt"] = {  "float","10","Timestep during random uplift forcing.","Used to perturb the landscape."};

  int_default_map["random_cycles"] = 4;
  help_map["random_cycles"] = {  "int","3","Number of different uplift rates in the random uplift forcing.","Used to perturb the landscape."};

  // some parameters for a spatially varying K and U
  bool_default_map["make_spatially_varying_K"] = false;   // This just prints a raster
  help_map["make_spatially_varying_K"] = {  "bool","false","This uses rastermaker to create a K raster and prints it","For running models with spatially variable K without a lithocube"};

  float_default_map["spatially_varying_max_K"] = 0.0001;
  help_map["spatially_varying_max_K"] = {  "float","0.0001","For making a random K field this is the maximum K value.","For running models with spatially variable K without a lithocube."};

  float_default_map["spatially_varying_min_K"] = 0.000001;
  help_map["spatially_varying_min_K"] = {  "float","0.000001","For making a random K field this is the minimum K value.","For running models with spatially variable K without a lithocube."};

  int_default_map["min_blob_size"] = 50;
  help_map["spatially_varying_min_K"] = {  "int","50","For making a random K field this is the minimum width of a random K blob in pixels (they are square).","For running models with spatially variable K without a lithocube."};

  int_default_map["max_blob_size"] = 100;
  help_map["max_blob_size"] = {  "int","100","For making a random K field this is the mmaximum width of a random K blob in pixels (they are square).","For running models with spatially variable K without a lithocube."};

  int_default_map["n_blobs"] = 10;
  help_map["n_blobs"] = {  "int","10","For making a random K field this is the number of blobs in the raster.","For running models with spatially variable K without a lithocube."};

  bool_default_map["spatially_varying_forcing"] = false;
  help_map["spatially_varying_forcing"] = {  "bool","false","Turn on one of spatially varying U or K","For running models with spatially variable U or K without a lithocube"};

  bool_default_map["spatially_varying_K"] = false;
  help_map["spatially_varying_K"] = {  "bool","false","The tells the model to read a K raster","For running models with spatially variable K without a lithocube"};

  bool_default_map["spatially_varying_U"] = false;
  help_map["spatially_varying_U"] = {  "bool","false","The tells the model to read a U raster","For running models with spatially variable U"};

  bool_default_map["calculate_K_from_relief"] = false;
  help_map["calculate_K_from_relief"] = {  "bool","false","This takes a DEM and an uplift rate and back calculates the K value for each pixel","For constraining the K values in a real DEM"};

  int_default_map["spatial_K_method"] = 0;
  help_map["spatial_K_method"] = {  "int","0","Method for making spatially variable K options are 0 for blobs 1 dinsnae work and 2 has half the raster at higher K","For variable K values"};

  int_default_map["spatial_U_method"] = 1;          // 0 doesn't work
  help_map["spatial_U_method"] = {  "int","1","Method for making spatially variable U options are 0 dinsnae work 1 sine wave and 2 has half the raster at higher U","For variable U values"};

  bool_default_map["load_K_raster"] = false;
  help_map["load_K_raster"] = {  "bool","false","Loads the K raster from file","For imposing spatially variable K"};

  bool_default_map["load_U_raster"] = false;
  help_map["load_U_raster"] = {  "bool","false","Loads the U raster from file","For imposing spatially variable U"};

  float_default_map["spatial_K_factor"] = 3;
  help_map["spatial_K_factor"] = {  "float","3","For K variation this the the factor over which K varies","For variable K values"};

  float_default_map["spatial_variation_time"] = 20000;
  help_map["spatial_variation_time"] = {  "float","20000","I don't understand what this is for and I wrote it--SMM","For variable K values"};

  bool_default_map["make_spatially_varying_U"] = false;  // This just prints a raster
  help_map["make_spatially_varying_U"] = {  "bool","false","This uses rastermaker to create a U raster and prints it","For running models with spatially variable U"};

  float_default_map["min_U_for_spatial_var"] = 0.0001;
  help_map["min_U_for_spatial_var"] = {  "float","0.0001","For spatially varied U this is the minimum value.","For running models with spatially variable U."};

  float_default_map["max_U_for_spatial_var"] = 0.0005;
  help_map["max_U_for_spatial_var"] = {  "float","0.0005","For spatially varied U this is the maximum value.","For running models with spatially variable U."};

  int_default_map["K_smoothing_steps"] = 2;
  help_map["K_smoothing_steps"] = {  "int","2","For spatially varied K this diffuses K around the edges of the blobs.","For running models with spatially variable K. ensures you don't get discontinuities at the boundaries."};

  float_default_map["spatial_dt"] = 100;
  help_map["spatial_dt"] = {  "float","100","Timestep when running spatially varied K or U.","For running models with spatially variable K. ensures you don't get discontinuities at the boundaries."};

  int_default_map["spatial_cycles"] = 5;
  help_map["spatial_cycles"] = {  "int","5","Number of cycles for varied K or U.","For running models with spatially variable K. Why did I do this and not just set the end time?--SMM."};

  bool_default_map["use_adaptive_timestep"] = true;
  help_map["use_adaptive_timestep"] = {  "bool","true","This enables the variable timestep.","With this enabled you can speed up the model or avoid numerical instability but in some cases it really slows it down."};

  float_default_map["maximum_timestep"] = 500;
  help_map["maximum_timestep"] = {  "float","500","For adaptive timestepping this caps the timestep.","Only used with use_adaptive_timestep."};

  bool_default_map["snap_to_steep_for_spatial_uplift"] = false;
  help_map["snap_to_steep_for_spatial_uplift"] = {  "bool","false","If you are using variable upflift this will snap either to minimum or maximum uplift.","Whether min or max depends on nap_to_minimum_uplift."};

  bool_default_map["snap_to_minimum_uplift"] = true;
  help_map["snap_to_minimum_uplift"] = {  "bool","true","If you are using variable upflift this will snap to minimum uplift but if false will snap to maximum uplift.","snap_to_steep_for_spatial_uplift needs to be true."};

  int_default_map["increase_amt"] = 10;
  help_map["increase_amt"] = {  "int","10","For variable K or U that have half the forcing field at a different value this is the multipler for the high value.","For say 10 the north will have 10x the uplift of the south."};

  bool_default_map["finish_with_steady_forcing"] = false;
  help_map["finish_with_steady_forcing"] = {  "bool","false","If you are using variable upflift this will finish the run with a period of steady uplift.","To smooth out the raster at the end. Although I never really used this--SMM"};



  // Some faults or base level falls
  bool_default_map["simple_fault_uplift"] = false;
  help_map["simple_fault_uplift"] = {  "bool","false","This just raises the north half of the raster by a fixed amount","Some basic transient forcing for toy models."};

  int_default_map["uplift_amt"] = 10;
  help_map["simple_fault_uplift"] = {  "int","10","For simple fault uplift this is the instantaneous uplift of the fault","Some basic transient forcing for toy models."};

  int_default_map["uplift_time"] = 50000;
  help_map["uplift_time"] = {  "int","50000","For simple fault uplift this is the time when the fault goes","Some basic transient forcing for toy models."};

  float_default_map["uplift_rate"] = 0.0002;
  help_map["uplift_rate"] = {  "float","0.0002","Uplift rate for various transient scenarios and snapping","Some basic transient forcing for toy models."};

  bool_default_map["base_level_fall"] = false;
  help_map["base_level_fall"] = {  "bool","false","This just raises drops one side of the raster by a fixed amoount","Some basic transient forcing for toy models."};

  int_default_map["transient_print_interval"] = 2;
  help_map["transient_print_interval"] = {  "int","2","Print interval for transient runs","This is controled by dt."};


  // An imposed transient scenario
  bool_default_map["run_transient_base_level"] = false;
  help_map["run_transient_base_level"] = {  "bool","true","Tests the transient component","Still under construction."};

  // These are paramaters for the transient base level file
  string_default_map["transient_channel_csv_name"] = "NULL";  // This contains the fixed channel data
  help_map["transient_channel_csv_name"] = {  "string","NULL","This is the name of a csv file with a fixed channel.","Channel needs column headers latitude longitude elevation and has csv in the name. Also has timesteps in columns"};

  bool_default_map["test_transient_channel"] = false;
  help_map["test_transient_channel"] = {  "bool","false","This checks the format of the transient channel and exits.","For bug checking."};

  bool_default_map["use_transient_outlet"] = false;
  help_map["use_transient_outlet"] = {  "bool","false","If true you read the transient_outlet_fname and and drive the transient model with baselevel fall.","This runs the transient model with outlet elevations only."};

  string_default_map["transient_outlet_fname"]  = "outlet_elevations.csv";
  help_map["transient_outlet_fname"] = {  "string","outlet_elevations.csv","This contains the filename for the outlet elevations through time if you are running the model with the outlet only.","Only used if use_transient_outlet is true."};

  string_default_map["baselevel_switch_column"] = "MEAS_Aare";
  help_map["baselevel_switch_column"] = {  "string","MEAS_Aare","This is the name of a column used to popuate the baselevel switch. If this colum has a zero the baselevel switch is 0. If it is not zero the switch is 1.","If the switch is 1 it means that the pixel will have a forced drainage area but elevation will evolve transiently. A 0 means it is a baselevel node with a fixed elevation."};

  string_default_map["transient_channel_timing_prefix"] = "t";  // This contains the fixed channel data
  help_map["transient_channel_timing_prefix"] = {  "string","t","The prefix of the timing column headers.","If timing columns are t5 t10 etc then the prefix would be t"};

  float_default_map["transient_channel_timing_multiplier"] = 1000;  // This contains the fixed channel data
  help_map["transient_channel_timing_multiplier"] = {  "float","1000","How many years each column integer represents.","I t5 is 5000 years then the multiplyer is 1000"};

  int_default_map["transient_channel_n_time_columns"] = 10;  // This contains the fixed channel data
  help_map["transient_channel_n_time_columns"] = {  "int","10","Number of transient time steps you will attempt to read.","This looks for column names so should be resilient to crashing."};

  int_default_map["transient_channel_phase_steps"] = 5;  // This contains the fixed channel data
  help_map["transient_channel_phase_steps"] = {  "int","5","How many timing multipliers each phase step you will attempt.","If timing columns are t5 t10 etc then the steps will be 5."};

  float_default_map["transient_maximum_time"] = 10500;
  help_map["transient_maximum_time"] = {  "float","10500","Sets the maximum time for the transient run.","To get the full baselevel record set this longer than the baselelvel steps"};

  bool_default_map["use_lithocube_for_transient"] = false;
  help_map["use_lithocube_for_transient"] = {  "bool","false","Uses the lithocube in the transient runs if true.","You need to vo file for this to work."};

  // tilt scenarios
  bool_default_map["run_instantaneous_tilt"] = false;
  help_map["run_instantaneous_tilt"] = {  "bool","false","A tilting scenario where you tilt the whole raster instantaneously.","Set tilt with tilt_angle and tilt_boundary."};


  bool_default_map["run_progressive_tilt"] = false;
  help_map["run_progressive_tilt"] = {  "bool","false","A tilting scenario where you tilt the raster progressively.","Set tilt with tilt_angle and tilt_boundary."};

  int_default_map["tilt_time"] = 100000;
  help_map["tilt_time"] = {  "int","100000","For run_progressive_tilt scenario this is the time over which the tile will occur.","Set tilt with tilt_angle and tilt_boundary."};

  float_default_map["tilt_angle"] = 10;
  help_map["tilt_angle"] = {  "float","10","For tilt scenarios this is the final angle of the tilt.","Enable with one of run_progressive_tilt or run_instantaneous_tilt."};

  string_default_map["tilt_boundary"] = "S";
  help_map["tilt_boundary"] = {  "string","S","For tilt scenarios this is the side of the raster that will serve as the hinge.","Enable with one of run_progressive_tilt or run_instantaneous_tilt."};

  int_default_map["tilt_print_interval"] = 150;
  help_map["tilt_print_interval"] = {  "int","150","For tilt scenarios this is the printing interval.","Enable with one of run_progressive_tilt or run_instantaneous_tilt."};
 

  // Various options to load and process lithologic data
  bool_default_map["only_test_lithology"] = false;
  help_map["only_test_lithology"] = {  "bool","false","Loads the lithology data and stops there if this is true.","Used to check lithocube files."};

  string_default_map["vo_filename"] = "test.vo";
  help_map["vo_filename"] = {  "string","test.vo","This is the file that contains metadata for the lithocube.","Lithocubes are made of a small .vo file and a massive .asc file"};

  string_default_map["vo_asc_filename"] = "test.asc";
  help_map["vo_asc_filename"] = {  "string","test.asc","This is the name of the lithocube file that has the actual three dimensional lithology types in it.","You need to also give it the .vo file that contains the metadata about this file."};

  bool_default_map["print_K_raster"] = false;
  help_map["print_K_raster"] = {  "bool","false","Prints the K values in the lithocube testing.","Generates a raster with the lithocube derived K values."};

  bool_default_map["print_lithocode_raster"] = false;
  help_map["print_lithocode_raster"] = {  "bool","false","Prints the lithocodes in the lithocube testing.","Generates a raster with the lithocube derived lithocodes."};

  bool_default_map["print_exhumation_and_cumulative_uplift"] = false;
  help_map["print_exhumation_and_cumulative_uplift"] = {  "bool","false","Prints the cumulative uplift and exhumation rasters for use with transient runs and the lithocube.","Generates rasters for checking exhumation and cumulative uplift."};


  // options for loading lookup tables
  bool_default_map["load_lithocode_to_K_csv"] = false;
  help_map["load_lithocode_to_K_csv"] = {  "bool","false","Reads the lookup table converting lithocode to K if true.","Filename designated by lookup_filename_K."};

  bool_default_map["load_lithocode_to_Sc_csv"] = false;
  help_map["load_lithocode_to_Sc_csv"] = {  "bool","false","Reads the lookup table converting lithocode to Sc if true.","Filename designated by lookup_filename_Sc."};

  string_default_map["lookup_filename_K"] = "lookup_table_strati_K.csv";
  help_map["lookup_filename_K"] = {  "string","lookup_table_strati_K.csv","The name of the lithocode to K csv file.","Column names are Strati and K."};

  string_default_map["lookup_filename_Sc"] = "lookup_table_strati_Sc.csv";
  help_map["lookup_filename_Sc"] = {  "string","lookup_table_strati_Sc.csv","The name of the lithocode to Sc csv file.","Column names are Strati and Sc."};


  // option to have forbidden lithocodes
  string_default_map["forbidden_lithocodes"] = "NULL";
  help_map["forbidden_lithocodes"] = {  "string","NULL","If a lithocube value takes one of these then it is set to a default value.","Used to override lithologic unit K and Sc values with a default."};
  
  //Option to adjust elevation of raster before lithocube snapping
  float_default_map["elevation_change"] = 0;
  help_map["elevation_change"] = {  "float","0","Adjusts the elevation of the base raster by this amount.","For when you need to lift or drop the base raster."};
  
  //=========================================================================
  //
  //.#####....####...#####....####...##...##..######..######..######..#####..
  //.##..##..##..##..##..##..##..##..###.###..##........##....##......##..##.
  //.#####...######..#####...######..##.#.##..####......##....####....#####..
  //.##......##..##..##..##..##..##..##...##..##........##....##......##..##.
  //.##......##..##..##..##..##..##..##...##..######....##....######..##..##.
  //
  //..####...##..##..######...####...##..##...####..                         
  //.##..##..##..##..##......##..##..##.##...##.....                         
  //.##......######..####....##......####.....####..                         
  //.##..##..##..##..##......##..##..##.##.......##.                         
  //..####...##..##..######...####...##..##...####..                         
  //============================================================================

  // Use the parameter parser to get the maps of the parameters required for the analysis
  LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
  map<string,float> this_float_map = LSDPP.get_float_parameters();
  map<string,int> this_int_map = LSDPP.get_int_parameters();
  map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
  map<string,string> this_string_map = LSDPP.get_string_parameters();

  if(f_name == "cry_for_help.txt")
  {
    cout << "I am going to print the help and exit." << endl;
    cout << "You can find the help in the file:" << endl;
    cout << "./MuddPILEDriver-README.csv" << endl;
    string help_prefix = "MuddPILEDriver-README";
    LSDPP.print_help(help_map, help_prefix, version_number, citation);
    exit(0);
  }


   
  list<int> forbidden_lithocodes;
  if (this_string_map["forbidden_lithocodes"] != "NULL")
  {
    string this_key = "forbidden_lithocodes";
    vector<int> temp_forbidden_codes;
    temp_forbidden_codes = LSDPP.parse_int_vector(this_key);

    cout << "Reading forbidden lithocodes. They are: " << endl;
    for (size_t i = 0; i<temp_forbidden_codes.size(); i++)
    {
      forbidden_lithocodes.push_back(temp_forbidden_codes[i]);
      cout << temp_forbidden_codes[i] << endl;
    }
  }


  // Now print the parameters for bug checking
  //LSDPP.print_parameters();

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
  //mod.print_parameters();

  // need this to keep track of the end time
  float current_end_time = 0;

  //============================================================================
  //.##......######..######..##..##...####...........##.......####....####...#####..
  //.##........##......##....##..##..##..##..........##......##..##..##..##..##..##.
  //.##........##......##....######..##..##..........##......##..##..######..##..##.
  //.##........##......##....##..##..##..##..........##......##..##..##..##..##..##.
  //.######..######....##....##..##...####...........######...####...##..##..#####..
  //============================================================================
  // Logic for loading lithology
  //============================================================================ 
  if(this_bool_map["only_test_lithology"])
  {
    cout << endl << endl << endl << "============================================" << endl;
    cout << "Let's load the lithology! I am doing this as a test so this routine" << endl;
    cout << "Will exit once the lithology is loaded." << endl;
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

    
    // Loading lithocode
    if( this_bool_map["load_lithocode_to_K_csv"])
    {   
      string this_lookup_table_K = DATA_DIR+this_string_map["lookup_filename_K"];
      cout << "I am going to load a lookup table for K" << endl;
      cout << "The filename is: " << this_lookup_table_K << endl;
      LSDLC.read_strati_to_K_csv(this_lookup_table_K);
    }
    else
    {
      LSDLC.hardcode_fill_strati_to_K_map();  
    }  
    // Let's create the map of stratigraphy codes and K values    
    if( this_bool_map["load_lithocode_to_Sc_csv"])
    {   
      string this_lookup_table_Sc = DATA_DIR+this_string_map["lookup_filename_Sc"];
      cout << "I am going to load a lookup table for Sc" << endl;
      cout << "The filename is: " << this_lookup_table_Sc << endl;
      LSDLC.read_strati_to_Sc_csv(this_lookup_table_Sc);
    }
    else
    {
      LSDLC.hardcode_fill_strati_to_Sc_map();  
    }  

    // // test printing of a layer
    // string test_layer_name = DATA_DIR+DEM_ID+"_Layer61";
    // int test_layer = 61;
    // LSDIndexRaster test_layer_raster = LSDLC.get_raster_from_layer(test_layer);
    // LSDRaster K_raster_from_lithocube_layer = LSDLC.index_to_K(test_layer_raster, 0.000001);
    // cout << "Printing a K raster from layer " << test_layer << endl;
    // K_raster_from_lithocube_layer.write_raster(test_layer_name, "bil");

    string lithocodes_raster_name = DATA_DIR+OUT_ID+"_LithoCodes";
    string test_K_raster_name = DATA_DIR+OUT_ID+"_KRasterLSDLC";
    string test_Sc_raster_name = DATA_DIR+OUT_ID+"_ScRasterLSDLC";

    // LSDLC.write_litho_raster_from_elevations(elev_raster, test_litho_raster_name, "bil");
    cout << "Getting the lithology codes" << endl;
    LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster,forbidden_lithocodes);
    LSDRaster lithocodes_raster(lithocodes_index_raster);
    cout << "Printing a lithocodes raster" << endl;
    lithocodes_raster.write_raster(lithocodes_raster_name, "bil");    
    cout << "Converting the litho codes to K values" << endl;
    LSDRaster K_raster_from_lithocube = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);
    cout << "Converting the litho codes to Sc values" << endl;
    LSDRaster Sc_raster_from_lithocube = LSDLC.index_to_Sc(lithocodes_index_raster, 0.21);

    cout << "Printing a K raster" << endl;
    K_raster_from_lithocube.write_raster(test_K_raster_name, "bil");
    cout << "Wrote a K raster for you" << endl;

    cout << "Printing an Sc raster" << endl;
    Sc_raster_from_lithocube.write_raster(test_Sc_raster_name, "bil");
    cout << "Wrote an Sc raster for you" << endl;

    //LSDLC.get_K_raster_from_raster(elev_raster);



    cout << "I've tested the lithology object and now am exiting." << endl;
    exit(0);
    

  } 

  //============================================================================
  //.######..##..##..######..######..######...####...##.....
  //...##....###.##....##......##......##....##..##..##.....
  //...##....##.###....##......##......##....######..##.....
  //...##....##..##....##......##......##....##..##..##.....
  //.######..##..##..######....##....######..##..##..######.
  //........................................................
  //.#####....####....####...######..######..#####..        
  //.##..##..##..##..##........##....##......##..##.        
  //.#####...######...####.....##....####....#####..        
  //.##..##..##..##......##....##....##......##..##.        
  //.##..##..##..##...####.....##....######..##..##. 
  //============================================================================
  // Logic for reading an initial surface
  // It will look for a raster but if it doesn't find one it will default back to
  // fixed rows and columns. If you do load a raster it will override any
  // instructions about NRows, NCols, and DataResolution
  //============================================================================
  LSDRaster InitialRaster;
  LSDRaster filled_topography;
  bool create_initial_surface = false;
  if(this_bool_map["read_initial_raster"])
  {
    cout << endl << endl << endl << "====================================" << endl;
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
      InitialRaster = temp_raster;

      map<string,string> GRS = temp_raster.get_GeoReferencingStrings();
      string CS = GRS["ENVI_coordinate_system"];
      string MI = GRS["ENVI_map_info"];
      cout << "I am loading a raster with Cooridinate string: " << endl;
      cout << CS << endl;
      cout << "and MI is: " << endl;
      cout << MI << endl;

      this_int_map["NRows"] = temp_raster.get_NRows();
      this_int_map["NCols"] = temp_raster.get_NCols();
      this_float_map["DataResolution"] = temp_raster.get_DataResolution();


      //==========================================================================
      // Fill the raster
      //==========================================================================
      LSDRaster carved_topography;  // filled topography is declared outside this logic
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


      // some tests for the single channel
      if( this_bool_map["impose_single_channel_printing"] || 
          this_bool_map["only_test_single_channel"] || 
          this_string_map["fixed_channel_csv_name"] != "NULL")
      { 
        cout << "You seem to have some logic for a single channel in this model run." << endl;
        cout << "I need to check if your flow distance and elevation columns are in the file " << endl;
        string fd_column_name = this_string_map["single_channel_fd_string"];
        string elevation_column_name = this_string_map["single_channel_elev_string"];

        LSDRasterInfo RI(filled_topography);
        // Get the latitude and longitude
        cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
        LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );
        single_channel_data.print_data_map_keys_to_screen();  

        if (not single_channel_data.is_column_in_csv(elevation_column_name))
        {
          string elevation_fragment = "elev";
          elevation_column_name = single_channel_data.find_string_in_column_name(elevation_fragment);
          if (elevation_column_name == "NULL")
          {
            cout << "Fatal error. I cannot find the elevation column name. " << endl;
            cout << "Check your single channel data file and the single_channel_elev_string. " << endl;
            exit(0);
          }

        }
        if (not single_channel_data.is_column_in_csv(fd_column_name))
        {
          string elevation_fragment = "flow";
          elevation_column_name = single_channel_data.find_string_in_column_name(elevation_fragment);
          if (elevation_column_name == "NULL")
          {
            cout << "Fatal error. I cannot find the flow distance column name. " << endl;
            cout << "Check your single channel data file and the single_channel_fd_string. " << endl;
            exit(0);
          }
        }

        cout << "You flow distance column name is " << fd_column_name << " and elev is: " << elevation_column_name << endl;
        cout << "If these are not in the column names this routine will fail!" << endl;
        cout << "WARNING: whitespace is removed from column names!!" << endl;
      }

      // .#####...##..##..######..######..######..#####..        
      // .##..##..##..##..##......##......##......##..##.        
      // .#####...##..##..####....####....####....#####..        
      // .##..##..##..##..##......##......##......##..##.        
      // .#####....####...##......##......######..##..##.        
      // ................................................        
      // ..####...##..##...####...##..##..##..##..######..##.....
      // .##..##..##..##..##..##..###.##..###.##..##......##.....
      // .##......######..######..##.###..##.###..####....##.....
      // .##..##..##..##..##..##..##..##..##..##..##......##.....
      // ..####...##..##..##..##..##..##..##..##..######..######.
      // This buffers and implements the single channel
      if (this_bool_map["buffer_single_channel"] && this_string_map["fixed_channel_csv_name"] != "NULL")
      {
        cout << endl << endl << endl <<"=================================================" << endl;
        cout << "Let me buffer the single channel." << endl;
        cout << "I am doing this to stop the channel from segmenting after snapping." << endl;
        cout << "Don't worry, it is only getting buffered by one pixel." << endl;
        string fd_column_name = this_string_map["single_channel_fd_string"];
        string elevation_column_name = this_string_map["single_channel_elev_string"];

        LSDRasterInfo RI(filled_topography);
        // Get the latitude and longitude
        cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
        LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );
        single_channel_data.print_data_map_keys_to_screen();

        cout << "Subtracting " << this_float_map["single_channel_drop"] << " m elevation from the single channel" << endl;
        single_channel_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
        single_channel_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
        // single_channel_data.print_data_to_csv(this_string_map["test_single_channel_name"]);

        LSDRasterMaker RM(filled_topography);

        if(this_bool_map["lower_raster_before_channel_buffer"])
        {
          // We subtract some elevation because the imposed channel has been dropped by 100m prior to saving as csv
          RM.add_value(-this_float_map["lower_raster_before_buffer_amount"]);
        }

        RM.impose_channels_with_buffer(single_channel_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);

        filled_topography = RM.return_as_raster();

        filled_topography.write_raster(OUT_DIR+OUT_ID+"_testbuffer","bil");

        // This is new as of 09/02/2021 and makes the model raster (for snapping etc) include the buffer
        temp_raster = filled_topography;
        cout << "========================================" << endl << endl << endl;
      }

      // This is a testing routine for the single channel
      if(this_bool_map["only_test_single_channel"])
      {
        cout << "Let me test the single channel modification routines." << endl;
        cout << "I am going to check the single channel and exit." << endl;
        cout << "If you want me to continue you need to set only_test_single_channel: false" << endl;
        cout << endl << endl << endl <<"=================================================" << endl;
        string fd_column_name = this_string_map["single_channel_fd_string"];
        string elevation_column_name = this_string_map["single_channel_elev_string"];

        LSDRasterInfo RI(filled_topography);
        // Get the latitude and longitude
        cout << "I am reading points from the file: "+ this_string_map["fixed_channel_csv_name"] << endl;
        LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );
        single_channel_data.print_data_map_keys_to_screen();

        cout << "Subtracting " << this_float_map["single_channel_drop"] << " m elevation from the single channel" << endl;
        single_channel_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
        single_channel_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
        //single_channel_data.print_data_to_csv(this_string_map["test_single_channel_name"]);


        single_channel_data.print_data_map_keys_to_screen();

        string initial_elev_column = "initial_elevation(m)";
        single_channel_data.burn_raster_data_to_csv(InitialRaster,initial_elev_column);

        single_channel_data.print_data_to_csv(this_string_map["test_single_channel_name"]);

        single_channel_data.print_row_and_col_to_csv("check_r_and_c.csv");

        exit(0);
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
        cout << "Let me extract the single channel for you." << endl;
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
      //mod.initialise_taper_edges_and_raise_raster(1);
      create_initial_surface = false;
      this_bool_map["force_dissect"] = false;


      if( this_bool_map["test_single_basin_outlet"])
      {
        cout << "I am going to test routing to a single outlet." << endl;
        LSDRasterInfo RI(filled_topography);
        LSDRasterMaker Test_basin_1(temp_raster);
        LSDRasterMaker Test_basin_2(temp_raster);
        LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        cout << "Buffering the long way" << endl;
        Test_basin_1.buffer_basin_to_single_outlet(single_channel_data, this_float_map["min_slope_for_fill"]);

        cout << "Buffering the short way" << endl;
        Test_basin_2.buffer_basin_to_single_outlet(this_float_map["min_slope_for_fill"]);

        LSDRaster Test1 = Test_basin_1.return_as_raster();
        //bool belowthresholdisnodata = true;
        //Test1 = Test1.mask_to_nodata_using_threshold(-100,belowthresholdisnodata);
        LSDRaster Test2 = Test_basin_2.return_as_raster();
        //Test2 = Test2.mask_to_nodata_using_threshold(-100,belowthresholdisnodata);

        string T1_fname = OUT_DIR+OUT_ID+"_basintest1";
        Test1.write_raster(T1_fname,raster_ext);

        string T2_fname = OUT_DIR+OUT_ID+"_basintest2";
        Test2.write_raster(T2_fname,raster_ext);

        cout << "Done with the test! Exiting" << endl;
        exit(0);

      }


      // .#####...######...####...##...##...####...##..##..#####..
      // .##..##....##....##..##..###.###..##..##..###.##..##..##.
      // .##..##....##....######..##.#.##..##..##..##.###..##..##.
      // .##..##....##....##..##..##...##..##..##..##..##..##..##.
      // .#####...######..##..##..##...##...####...##..##..#####..
      // .........................................................
      // ..####....####...##..##...####...#####...######......... 
      // .##......##..##..##..##..##..##..##..##..##............. 
      // ..####...##.###..##..##..######..#####...####........... 
      // .....##..##..##..##..##..##..##..##..##..##............. 
      // ..####....#####...####...##..##..##..##..######......... 
      // ........................................................ 
      // .#####...######..#####...##.......####....####...######. 
      // .##..##..##......##..##..##......##..##..##..##..##..... 
      // .#####...####....#####...##......######..##......####... 
      // .##..##..##......##......##......##..##..##..##..##..... 
      // .##..##..######..##......######..##..##...####...######. 
      // ........................................................ 
      // now for logic if you want a diamond square replacement
      if (this_bool_map["convert_initial_raster_to_diamond_square_fractal"])
      {
        int this_frame;
        cout << endl << endl << endl << "=<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]<>[]===" << endl;
        cout << "I going to replace your raster with a diamond square fractal surface." << endl;
        cout <<"then dissect the landscape for you." << endl;

        int smallest_dimension;
        if( mod.get_NRows() <= mod.get_NCols())
        {
          smallest_dimension = mod.get_NRows();
        }
        else
        {
          smallest_dimension = mod.get_NCols();
        }

        cout << "Getting the diamond square dimensions for your loaded DEM. " << endl;
        cout << "The DEM has " << mod.get_NRows() << " rows and " << mod.get_NCols() << " columns." << endl;
        // now get the largest possible dimension
        float lps = floor( log2( float(smallest_dimension) ) );
        int largest_possible_scale = int(lps);
        if( this_int_map["diamond_square_feature_order"] > largest_possible_scale )
        {
          this_int_map["diamond_square_feature_order"] = largest_possible_scale;
        }
        cout << "    The largest dimension is: " << smallest_dimension << " pixels." << endl;
        cout << "    So the biggest scale is: " << pow(2,largest_possible_scale) << " pixels or 2 to the " << largest_possible_scale << " power" << endl;
        cout << "    I am setting the scale to: " << this_int_map["diamond_square_feature_order"] << endl;

        // Now actually build the fractal surface
        mod.intialise_diamond_square_fractal_surface(this_int_map["diamond_square_feature_order"], this_float_map["diamond_square_relief"]);

        //cout << "   Let me taper the edges of this fractal surface for you so that everything drains to the edge." << endl;
        //cout << "   I am tapering along the " << this_int_map["taper_rows"] << " row closes to the N and S boundaries" << endl;
        //mod.initialise_taper_edges_and_raise_raster(this_int_map["taper_rows"]);

        //cout << "   I am superimposing a parabola with a relief of "  << this_float_map["parabola_relief"] << " metres" << endl;
        //mod.superimpose_parabolic_surface(this_float_map["parabola_relief"]);

        cout << "   I am going to fill and then roughen the surface for you. Roughness elements will " << endl;
        cout << "   have a maximum amplitude of " << this_float_map["roughness_relief"]<< " metres." << endl;
        mod.raise_and_fill_raster();
        mod.set_noise(this_float_map["roughness_relief"]);
        mod.random_surface_noise();
        mod.raise_and_fill_raster();

        this_frame = 9999;
        mod.print_rasters_and_csv( this_frame );

        // Now impose nodata
        cout << "Imposing nodata around the edges of your DS surface" << endl;
        mod.mask_to_nodata_impose_nodata(temp_raster);
        this_frame = 9998;
        mod.print_rasters_and_csv( this_frame );

        // Now tilt the edge up
        cout << "Tilting up the edge of your DS surface, and getting it to drain to a single outlet." << endl;
        LSDRasterInfo RI(filled_topography);
        LSDRaster this_raster = mod.return_as_raster();
        LSDRasterMaker new_basin(this_raster);
        LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        cout << "I am enforcing the slope on the single channel" << endl;
        string fd_column_name = this_string_map["single_channel_fd_string"];
        string elevation_column_name = this_string_map["single_channel_elev_string"];
        single_channel_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);

        // get the maximum elevation in this collection of points
        vector<float> elevs = single_channel_data.data_column_to_float(this_string_map["single_channel_elev_string"]);
        float max_elev = -99999;
        for (int i = 0; i< int(elevs.size()); i++)
        {
          if (elevs[i] > max_elev)
          {
            max_elev = elevs[i];
          }
        }
        cout << "The maximum elevation in the single channel is: " << max_elev << endl;

        new_basin.buffer_basin_to_single_outlet(single_channel_data, this_float_map["min_slope_for_fill"]);
        this_raster = new_basin.return_as_raster();
        mod.set_raster_data(this_raster);

        // uplift the raster to the maximum elevation
        cout << "I'm adding a fixed elevation to your DEM" << endl;
        mod.add_fixed_elevation(max_elev);

        this_frame = 9997;
        mod.print_rasters_and_csv( this_frame );

        // impose the channels then raise and fill
        mod.impose_channels(single_channel_data, this_string_map["single_channel_elev_string"]);
        mod.raise_and_fill_raster();

        this_frame = 9996;
        mod.print_rasters_and_csv( this_frame );

        mod.set_hillslope(false);
        //mod.set_timeStep( 1 );
        int this_uplift_mode = 0;     // block uplift
        float this_K = 0.001;
        float this_U = 0.0005;
        mod.set_K(this_K);

        // we run a snap just to see what happens
        cout << "Finished with the initial DS surface, snapping to 500m relief. " << endl;
        float desired_relief = 500;
        mod.fluvial_snap_to_steady_state_tune_K_for_relief(this_U, desired_relief);
        this_frame = 9995;
        mod.print_rasters_and_csv( this_frame );     

        // now do some repeat snapping
        int n_snaps = 20;
        long seed = time(NULL);
        mod.set_noise(2);
        for (int snap = 0; snap<n_snaps; snap++)
        {
          cout << "Snapping " << snap << " of " << n_snaps << endl;
          mod.random_surface_noise(seed);  
          mod.impose_channels(single_channel_data, this_string_map["single_channel_elev_string"]);
          mod.raise_and_fill_raster();  
          mod.fluvial_snap_to_steady_state_tune_K_for_relief(this_U, desired_relief);           
        }

        this_frame = 9994;
        mod.print_rasters_and_csv( this_frame );            


        // Now we do some fluvial incision
        cout << "Now I am going to run some fluvial incision at a small timestep" << endl;
        cout <<"I am going to use a moderate K value and uplift, since I don't" << endl;
        cout << "want to overexcavate the surface" << endl;
        mod.set_uplift( this_uplift_mode, this_U );
        this_frame = 1;
        int n_cycles = 5;
        for (int cycle = 1; cycle<=n_cycles; cycle++)
        {
          cout << "Cycle " << cycle << " of " << n_cycles << endl;
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
        this_frame = 9993;
        mod.print_rasters_and_csv( this_frame );

        // we run a snap just to see what happens
        cout << "Finished with the initial DS surface, snapping to 500m relief. " << endl;
        mod.raise_and_fill_raster();
        mod.fluvial_snap_to_steady_state_tune_K_for_relief(this_U, desired_relief);
        this_frame = 9992;
        mod.print_rasters_and_csv( this_frame );        


        current_end_time = 0;        
      } // Finished routine for replacing the current DEM with a diamond square 

      // This prints the channel network. Uses the chi tool so that we have the source pixels
      // As a column in the extraction
      if (this_bool_map["print_initial_channel_network"])
      {
        cout << "I am printing the initial channel network." << endl;
        if (this_bool_map["write_hillshade"])
        {
          cout << "Let me print the hillshade of the initial map for you. " << endl;
          float hs_azimuth = 315;
          float hs_altitude = 45;
          float hs_z_factor = 1;
          LSDRaster hs_raster = filled_topography.hillshade(hs_altitude,hs_azimuth,hs_z_factor);

          string hs_fname = OUT_DIR+OUT_ID+"_hs";
          hs_raster.write_raster(hs_fname,raster_ext);
        }
        
        cout << "\t Flow routing..." << endl;
        // get a flow info object
        LSDFlowInfo FlowInfo(boundary_conditions,filled_topography);

        // calculate the flow accumulation
        cout << "\t Calculating flow accumulation (in pixels)..." << endl;
        LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

        cout << "\t Converting to flow area..." << endl;
        LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

        if (this_bool_map["print_DrainageArea_raster"])
        {
          string DA_raster_name = OUT_DIR+OUT_ID+"_DArea";
          DrainageArea.write_raster(DA_raster_name,raster_ext);
        }

        // calculate the distance from outlet
        cout << "\t Calculating flow distance..." << endl;
        LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

        cout << "\t Loading Sources..." << endl;
        cout << "\t Source file is... " << CHeads_file << endl;

        // load the sources
        vector<int> sources;
        if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
        {
          cout << endl << endl << endl << "==================================" << endl;
          cout << "The channel head file is null. " << endl;
          cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
          sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

          cout << "The number of sources is: " << sources.size() << endl;
        }
        else
        {
          cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
          sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
          cout << "\t Got sources!" << endl;
        }

        // now get the junction network
        LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

        // get the chi coordinate
        float A_0 = 1.0;
        float thresh_area_for_chi = float( this_int_map["threshold_contributing_pixels"] ); // needs to be smaller than threshold channel area

        float movern = this_float_map["m"]/this_float_map["n"];
        LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

        vector<int> BaseLevelJunctions = JunctionNetwork.get_BaseLevelJunctions();

        cout << endl << endl << endl << "=================" << endl;
        for (int i = 0; i< int(BaseLevelJunctions.size()); i++)
        {
          cout << "bl["<<i<<"]: " << BaseLevelJunctions[i] << endl;
        }
        cout <<"=================" << endl;

        vector<int> source_nodes;
        vector<int> outlet_nodes;
        vector<int> baselevel_node_of_each_basin;
        int n_nodes_to_visit = 10;
        JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                    source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

        LSDChiTools ChiTool_chi_checker(FlowInfo);
        ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                            filled_topography, DistanceFromOutlet,
                            DrainageArea, chi_coordinate);

        string chi_data_maps_string = OUT_DIR+OUT_ID+"_initial_chi_data_map.csv";
        ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);

        if ( this_bool_map["convert_csv_to_geojson"])
        {
          string gjson_name = OUT_DIR+OUT_ID+"_initial_chi_data_map.geojson";
          LSDSpatialCSVReader thiscsv(chi_data_maps_string);
          thiscsv.print_data_to_geojson(gjson_name);
        }

      }   // Finished routine for printing the initial channel network

    }
    else
    {
      cout << "Header doesn't exist. I am switching to a fixed rows and columns." << endl;
      cout << "I will also create an initial surface for you." << endl;
      create_initial_surface = true;
    }


  }      // This is the end of the logic for reading the inital raster
  else
  {
    cout << "You have chosen not to read an initial raster so I will create an initial surface for you." << endl;
    create_initial_surface = true;
  }




  // a bit of logic to override spinup logic if you choose the diamond square spinup
  if(this_bool_map["diamond_square_spinup"])
  {
    cout << "You have chosen the diamond square spinup. This overrides all other spinup otions" << endl;
    cout << "It also overrides initial surface options" << endl;
    create_initial_surface = false;
    this_bool_map["cyclic_spinup"] = false;
    this_bool_map["spinup"] = false;
  }


  //.######..##..##..######..######..######...####...##.....        
  //...##....###.##....##......##......##....##..##..##.....        
  //...##....##.###....##......##......##....######..##.....        
  //...##....##..##....##......##......##....##..##..##.....        
  //.######..##..##..######....##....######..##..##..######.        
  //........................................................        
  //..####...##..##..#####...######...####....####...######.........
  //.##......##..##..##..##..##......##..##..##..##..##.............
  //..####...##..##..#####...####....######..##......####...........
  //.....##..##..##..##..##..##......##..##..##..##..##.............
  //..####....####...##..##..##......##..##...####...######.........
  //................................................................
  //..####...#####...######..##..##..##..##..#####..                
  //.##......##..##....##....###.##..##..##..##..##.                
  //..####...#####.....##....##.###..##..##..#####..                
  //.....##..##........##....##..##..##..##..##.....                
  //..####...##......######..##..##...####...##.....     
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

  if(this_bool_map["make_constant_Sc_raster"])
  {
    LSDRaster Temp_raster = mod.return_as_raster();
    cout << "I am going to make a contant raster for the Sc parameter." << endl;
    cout << "Setting K to "  << this_float_map["rudimentary_critical_slope"] << endl;
    LSDRasterMaker KRaster(Temp_raster);
    KRaster.set_to_constant_value(this_float_map["rudimentary_critical_slope"]);

    // write the raster
    string K_fname = OUT_DIR+OUT_ID+"_ConstScRaster";
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

    LSDRasterMaker URaster(Temp_raster);
    URaster.set_to_constant_value(this_float_map["rudimentary_steady_forcing_uplift"]);

    // write the raster
    string U_fname = OUT_DIR+OUT_ID+"_ConstURaster";
    string bil_name = "bil";

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



  //=============================================================================
  //
  //  /$$$$$$                                          /$$                    
  // /$$__  $$                                        |__/                    
  // | $$  \__/ /$$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$  /$$ /$$$$$$$   /$$$$$$ 
  // |  $$$$$$ | $$__  $$ |____  $$ /$$__  $$ /$$__  $$| $$| $$__  $$ /$$__  $$
  //  \____  $$| $$  \ $$  /$$$$$$$| $$  \ $$| $$  \ $$| $$| $$  \ $$| $$  \ $$
  //  /$$  \ $$| $$  | $$ /$$__  $$| $$  | $$| $$  | $$| $$| $$  | $$| $$  | $$
  // |  $$$$$$/| $$  | $$|  $$$$$$$| $$$$$$$/| $$$$$$$/| $$| $$  | $$|  $$$$$$$
  //  \______/ |__/  |__/ \_______/| $$____/ | $$____/ |__/|__/  |__/ \____  $$
  //                               | $$      | $$                     /$$  \ $$
  //                               | $$      | $$                    |  $$$$$$/
  //                               |__/      |__/                     \______/ 
  // 
  //==============================================================================                                    
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

    if (this_float_map["single_channel_drop"] != 0)
    {
      cout << endl << endl << endl << "==================================" << endl;
      cout << "I am dropping and adjusting the slope of the single channel." << endl;
      cout << "The data keys are: " << endl;
      source_points_data.print_data_map_keys_to_screen();
      cout << "==================================" << endl << endl << endl <<endl;

      string elevation_column_name = this_string_map["single_channel_elev_string"];
      string fd_column_name = this_string_map["single_channel_fd_string"];

      source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
      source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
    }
    else
    {
      string elevation_column_name = this_string_map["single_channel_elev_string"];
      string fd_column_name = this_string_map["single_channel_fd_string"];      
      source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
    }

    
    if(this_bool_map["print_updated_fixed_channel"])
    {
      cout << "I am going to print the updated fixed channel that has both an imposed slope and a possible excavation depth." << endl;
      source_points_data.print_data_to_csv(this_string_map["test_single_channel_name"]);
    }


    if (this_bool_map["only_test_impose_single_channel"])
    {
      cout << "=====================================" << endl;
      cout << "I'm not going to snap, I'm just imposing the channels" << endl;
      string csv_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.csv";
      source_points_data.print_data_to_csv(csv_name);
      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.geojson";
        source_points_data.print_data_to_geojson(gjson_name);
      }   

      LSDRaster testing_topo_raster = mod.return_as_raster();
      LSDRasterMaker RM(testing_topo_raster);
      RM.impose_channels_with_buffer(source_points_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);

      LSDRaster imposed_raster = RM.return_as_raster();
      string irastername = OUT_DIR+OUT_ID+"_imposed";
      imposed_raster.write_raster(irastername,"bil");

      exit(0);
    }

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
      cout << "The U value at (500,500) is: " << U_values.get_data_element(500,500) << endl;
      mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"], this_string_map["single_channel_elev_string"]);

      if (this_bool_map["use_initial_topography_as_max_elevation"])
      {
        mod.cap_elevations(InitialRaster);
      }

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
    // FINISHED SNAPPING


    if (this_bool_map["print_final_channel_network"])
    {
      LSDRaster ft,ct;


      if(this_bool_map["impose_single_channel_printing"])
      {
        // load the single channel
        LSDRasterInfo RI(InitialRaster);
        // Get the latitude and longitude
        //cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        if(this_bool_map["fixed_channel_dig_and_slope"])
        {
          string fd_column_name = this_string_map["single_channel_fd_string"];
          string elevation_column_name = this_string_map["single_channel_elev_string"];

          source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
          source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);  

          string csv_name = OUT_DIR+OUT_ID+"final_single_channel.csv";
          source_points_data.print_data_to_csv(csv_name);
          if ( this_bool_map["convert_csv_to_geojson"])
          {
            string gjson_name = OUT_DIR+OUT_ID+"final_single_channel.geojson";
            source_points_data.print_data_to_geojson(gjson_name);
          }            
        }

        mod.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);
      }

      LSDRaster new_muddpile_topo_raster = mod.return_as_raster(); 
      if(this_bool_map["carve_before_fill"])
      {
        ct = new_muddpile_topo_raster.Breaching_Lindsay2016();
        ft = ct.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        ft = new_muddpile_topo_raster.fill(this_float_map["min_slope_for_fill"]);
      }

      // Write the raster
      string R_name = OUT_DIR+OUT_ID+"Final_DEM_fill";
      ft.write_raster(R_name,"bil");
      R_name = OUT_DIR+OUT_ID+"Final_DEM";
      new_muddpile_topo_raster.write_raster(R_name,"bil");
  
      cout << "\t Flow routing..." << endl;
      LSDFlowInfo FlowInfo(boundary_conditions,ft);

      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      cout << "\t Converting to flow area..." << endl;
      LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

      // calculate the distance from outlet
      cout << "\t Calculating flow distance..." << endl;
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

      cout << "\t Loading Sources..." << endl;
      cout << "\t Source file is... " << CHeads_file << endl;

      // load the sources
      vector<int> sources;
      if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
      {
        cout << endl << endl << endl << "==================================" << endl;
        cout << "The channel head file is null. " << endl;
        cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
        sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

        cout << "The number of sources is: " << sources.size() << endl;
      }
      else
      {
        cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
        sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
        cout << "\t Got sources!" << endl;
      }

      // now get the junction network
      LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

      // get the chi coordinate
      float A_0 = 1.0;
      float thresh_area_for_chi = float( this_int_map["threshold_contributing_pixels"] ); // needs to be smaller than threshold channel area

      float movern = this_float_map["m"]/this_float_map["n"];
      LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

      vector<int> BaseLevelJunctions = JunctionNetwork.get_BaseLevelJunctions();

        cout << endl << endl << endl << "=================" << endl;
        for (int i = 0; i< int(BaseLevelJunctions.size()); i++)
        {
          cout << "bl["<<i<<"]: " << BaseLevelJunctions[i] << endl;
        }
        cout <<"=================" << endl;

      vector<int> source_nodes;
      vector<int> outlet_nodes;
      vector<int> baselevel_node_of_each_basin;
      int n_nodes_to_visit = 10;
      JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

      LSDChiTools ChiTool_chi_checker(FlowInfo);
      ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                          ft, DistanceFromOutlet,
                          DrainageArea, chi_coordinate);

      string chi_data_maps_string = OUT_DIR+OUT_ID+"_final_chi_data_map.csv";
      ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_final_chi_data_map.geojson";
        LSDSpatialCSVReader thiscsv(chi_data_maps_string);
        thiscsv.print_data_to_geojson(gjson_name);
      }  
    }

    cout << "Thanks for snapping!" << endl;
  }

  //=============================================================================
  //
  //  /$$$$$$                                          /$$                    
  // /$$__  $$                                        |__/                    
  // | $$  \__/ /$$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$  /$$ /$$$$$$$   /$$$$$$ 
  // |  $$$$$$ | $$__  $$ |____  $$ /$$__  $$ /$$__  $$| $$| $$__  $$ /$$__  $$
  //  \____  $$| $$  \ $$  /$$$$$$$| $$  \ $$| $$  \ $$| $$| $$  \ $$| $$  \ $$
  //  /$$  \ $$| $$  | $$ /$$__  $$| $$  | $$| $$  | $$| $$| $$  | $$| $$  | $$
  // |  $$$$$$/| $$  | $$|  $$$$$$$| $$$$$$$/| $$$$$$$/| $$| $$  | $$|  $$$$$$$
  //  \______/ |__/  |__/ \_______/| $$____/ | $$____/ |__/|__/  |__/ \____  $$
  //                               | $$      | $$                     /$$  \ $$
  //                               | $$      | $$                    |  $$$$$$/
  //                               |__/      |__/                     \______/ 
  // 
  //   /$$$$$$            /$$   /$$              /$$                              
  //  /$$__  $$          |__/  | $$             | $$                              
  // | $$  \__/  /$$$$$$  /$$ /$$$$$$   /$$$$$$$| $$  /$$$$$$   /$$$$$$   /$$$$$$ 
  // | $$       /$$__  $$| $$|_  $$_/  /$$_____/| $$ /$$__  $$ /$$__  $$ /$$__  $$
  // | $$      | $$  \__/| $$  | $$   |  $$$$$$ | $$| $$  \ $$| $$  \ $$| $$$$$$$$
  // | $$    $$| $$      | $$  | $$ /$$\____  $$| $$| $$  | $$| $$  | $$| $$_____/
  // |  $$$$$$/| $$      | $$  |  $$$$//$$$$$$$/| $$|  $$$$$$/| $$$$$$$/|  $$$$$$$
  //  \______/ |__/      |__/   \___/ |_______/ |__/ \______/ | $$____/  \_______/
  //                                                          | $$                
  //                                                          | $$                
  //                                                          |__/                
  // 
  //=============================================================================  
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


      if (this_float_map["single_channel_drop"] != 0)
      {
        cout << endl << endl << endl << "==================================" << endl;
        cout << "I am dropping and adjusting the slope of the single channel." << endl;
        cout << "The data keys are: " << endl;
        source_points_data.print_data_map_keys_to_screen();
        cout << "==================================" << endl << endl << endl <<endl;

        string elevation_column_name = this_string_map["single_channel_elev_string"];
        string fd_column_name = this_string_map["single_channel_fd_string"];

        source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
        source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
      }
      else
      {
        string elevation_column_name = this_string_map["single_channel_elev_string"];
        string fd_column_name = this_string_map["single_channel_fd_string"];      
        source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
      }


      if (this_bool_map["only_test_impose_single_channel"])
      {
        cout << "=====================================" << endl;
        cout << "I'm not going to snap, I'm just imposing the channels" << endl;
        string csv_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.csv";
        source_points_data.print_data_to_csv(csv_name);
        if ( this_bool_map["convert_csv_to_geojson"])
        {
          string gjson_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.geojson";
          source_points_data.print_data_to_geojson(gjson_name);
        }   

        LSDRaster testing_topo_raster = mod.return_as_raster();
        LSDRasterMaker RM(testing_topo_raster);
        RM.impose_channels_with_buffer(source_points_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);

        LSDRaster imposed_raster = RM.return_as_raster();
        string irastername = OUT_DIR+OUT_ID+"_imposed";
        imposed_raster.write_raster(irastername,"bil");

        exit(0);
      }


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

        if (this_bool_map["use_initial_topography_as_max_elevation"])
        {
          mod.cap_elevations(InitialRaster);
        }


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
      // FINISHED SNAPPING
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

      if (this_bool_map["use_initial_topography_as_max_elevation"])
      {
        mod.cap_elevations(InitialRaster);
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


    if (this_bool_map["print_final_channel_network"])
    {
      LSDRaster ft,ct;

      if(this_bool_map["impose_single_channel_printing"])
      {
        // load the single channel
        LSDRasterInfo RI(InitialRaster);
        // Get the latitude and longitude
        //cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        if(this_bool_map["fixed_channel_dig_and_slope"])
        {
          string fd_column_name = this_string_map["single_channel_fd_string"];
          string elevation_column_name = this_string_map["single_channel_elev_string"];

          source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
          source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);  

          string csv_name = OUT_DIR+OUT_ID+"final_single_channel.csv";
          source_points_data.print_data_to_csv(csv_name);
          if ( this_bool_map["convert_csv_to_geojson"])
          {
            string gjson_name = OUT_DIR+OUT_ID+"final_single_channel.geojson";
            source_points_data.print_data_to_geojson(gjson_name);
          }   
        }

        mod.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);
      }


      LSDRaster new_muddpile_topo_raster = mod.return_as_raster(); 
      if(this_bool_map["carve_before_fill"])
      {
        ct = new_muddpile_topo_raster.Breaching_Lindsay2016();
        ft = ct.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        ft = new_muddpile_topo_raster.fill(this_float_map["min_slope_for_fill"]);
      }

      // Write the raster
      string R_name = OUT_DIR+OUT_ID+"Final_DEM_fill";
      ft.write_raster(R_name,"bil");
      R_name = OUT_DIR+OUT_ID+"Final_DEM";
      new_muddpile_topo_raster.write_raster(R_name,"bil");
  
      cout << "\t Flow routing..." << endl;
      LSDFlowInfo FlowInfo(boundary_conditions,ft);

      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      cout << "\t Converting to flow area..." << endl;
      LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

      // calculate the distance from outlet
      cout << "\t Calculating flow distance..." << endl;
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

      cout << "\t Loading Sources..." << endl;
      cout << "\t Source file is... " << CHeads_file << endl;

      // load the sources
      vector<int> sources;
      if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
      {
        cout << endl << endl << endl << "==================================" << endl;
        cout << "The channel head file is null. " << endl;
        cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
        sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

        cout << "The number of sources is: " << sources.size() << endl;
      }
      else
      {
        cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
        sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
        cout << "\t Got sources!" << endl;
      }

      // now get the junction network
      LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

      // get the chi coordinate
      float A_0 = 1.0;
      float thresh_area_for_chi = float( this_int_map["threshold_contributing_pixels"] ); // needs to be smaller than threshold channel area

      float movern = this_float_map["m"]/this_float_map["n"];
      LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

      vector<int> BaseLevelJunctions = JunctionNetwork.get_BaseLevelJunctions();

      cout << endl << endl << endl << "=================" << endl;
      for (int i = 0; i< int(BaseLevelJunctions.size()); i++)
      {
        cout << "bl["<<i<<"]: " << BaseLevelJunctions[i] << endl;
      }
      cout <<"=================" << endl;

      vector<int> source_nodes;
      vector<int> outlet_nodes;
      vector<int> baselevel_node_of_each_basin;
      int n_nodes_to_visit = 10;
      JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

      LSDChiTools ChiTool_chi_checker(FlowInfo);
      ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                          ft, DistanceFromOutlet,
                          DrainageArea, chi_coordinate);

      string chi_data_maps_string = OUT_DIR+OUT_ID+"_final_chi_data_map.csv";
      ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_final_chi_data_map.geojson";
        LSDSpatialCSVReader thiscsv(chi_data_maps_string);
        thiscsv.print_data_to_geojson(gjson_name);
      }  
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






  //=======================================================================================
  //  /$$       /$$   /$$     /$$                                                        
  // | $$      |__/  | $$    | $$                                                        
  // | $$       /$$ /$$$$$$  | $$$$$$$   /$$$$$$   /$$$$$$$ /$$$$$$$   /$$$$$$   /$$$$$$ 
  // | $$      | $$|_  $$_/  | $$__  $$ /$$__  $$ /$$_____/| $$__  $$ |____  $$ /$$__  $$
  // | $$      | $$  | $$    | $$  \ $$| $$  \ $$|  $$$$$$ | $$  \ $$  /$$$$$$$| $$  \ $$
  // | $$      | $$  | $$ /$$| $$  | $$| $$  | $$ \____  $$| $$  | $$ /$$__  $$| $$  | $$
  // | $$$$$$$$| $$  |  $$$$/| $$  | $$|  $$$$$$/ /$$$$$$$/| $$  | $$|  $$$$$$$| $$$$$$$/
  // |________/|__/   \___/  |__/  |__/ \______/ |_______/ |__/  |__/ \_______/| $$____/ 
  //                                                                           | $$      
  //                                                                           | $$      
  //                                                                           |__/      
  //=======================================================================================
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
    if (this_float_map["elevation_change"] != 0)
    {
      elev_raster.AdjustElevation(this_float_map["elevation_change"]);
      cout << "Adjusted the elevation of your raster by " << this_float_map["elevation_change"] << endl;
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
    if( this_bool_map["load_lithocode_to_K_csv"])
    {   
      string this_lookup_table_K = DATA_DIR+this_string_map["lookup_filename_K"];
      cout << "I am going to load a lookup table for K" << endl;
      cout << "The filename is: " << this_lookup_table_K << endl;
      LSDLC.read_strati_to_K_csv(this_lookup_table_K);
    }
    else
    {
      LSDLC.hardcode_fill_strati_to_K_map();  
    }  
    


    // Now we'll get the lithology codes, convert them to K values and make a raster of K values for the model to use
    cout << "Getting the lithology codes" << endl;
    LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster,forbidden_lithocodes);
    cout << "Converting the litho codes to K values" << endl;
    LSDRaster K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);
       
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

    if (this_float_map["single_channel_drop"] != 0)
    {
      cout << endl << endl << endl << "==================================" << endl;
      cout << "I am dropping and adjusting the slope of the single channel." << endl;
      cout << "The data keys are: " << endl;
      source_points_data.print_data_map_keys_to_screen();
      cout << "==================================" << endl << endl << endl <<endl;

      string elevation_column_name = this_string_map["single_channel_elev_string"];
      string fd_column_name = this_string_map["single_channel_fd_string"];

      source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
      source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
    }
    else
    {
      string elevation_column_name = this_string_map["single_channel_elev_string"];
      string fd_column_name = this_string_map["single_channel_fd_string"];      
      source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
    }

    if (this_bool_map["only_test_impose_single_channel"])
    {
      cout << "=====================================" << endl;
      cout << "I'm not going to snap, I'm just imposing the channels" << endl;
      string csv_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.csv";
      source_points_data.print_data_to_csv(csv_name);
      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.geojson";
        source_points_data.print_data_to_geojson(gjson_name);
      }   

      LSDRaster testing_topo_raster = mod.return_as_raster();
      LSDRasterMaker RM(testing_topo_raster);
      RM.impose_channels_with_buffer(source_points_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);

      LSDRaster imposed_raster = RM.return_as_raster();
      string irastername = OUT_DIR+OUT_ID+"_imposed_buffer";
      imposed_raster.write_raster(irastername,"bil");

      exit(0);
    }

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
      cout << "The U value at (500,500) is: " << U_values.get_data_element(500,500) << endl;

      cout << "*****************************************" << endl;
      cout << "This is snap " << i+1 << " of " << n_snaps << endl;
      cout << "*****************************************" << endl;
      
      LSDRaster lithocodes_raster;

      // if(i>0)
      // {
      cout << "Let me create a raster of the topography from the previous loop." << endl;
      LSDRaster new_muddpile_topo_raster = mod.return_as_raster();
      cout << "Now I'll check where in the lithocube we are and make a raster of lithocodes." << endl;
      LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(new_muddpile_topo_raster,forbidden_lithocodes);
      cout << "I am now converting the lithocodes index raster to an LSDRaster" << endl;
      lithocodes_raster = LSDRaster(lithocodes_index_raster);
      cout << "Converting the litho codes to K values" << endl;
      K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);

      // }

      mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"], this_string_map["single_channel_elev_string"]);


      if (this_bool_map["use_initial_topography_as_max_elevation"])
      {
        mod.cap_elevations(InitialRaster);
      }

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
    // FINISHED SNAPPING

    if (this_bool_map["print_final_channel_network"])
    {
      cout << endl << endl << endl << "================" << endl; 
      cout << "We have got to the end of a lithosnappung run!" << endl;
      cout << "I am now going to print the final channel network" << endl;


      LSDRaster ft,ct;

      if(this_bool_map["impose_single_channel_printing"])
      {
        // load the single channel
        LSDRasterInfo RI(InitialRaster);
        // Get the latitude and longitude
        //cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        if(this_bool_map["fixed_channel_dig_and_slope"])
        {
          string fd_column_name = this_string_map["single_channel_fd_string"];
          string elevation_column_name = this_string_map["single_channel_elev_string"];

          source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
          source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);  

          string csv_name = OUT_DIR+OUT_ID+"final_single_channel.csv";
          source_points_data.print_data_to_csv(csv_name);
          if ( this_bool_map["convert_csv_to_geojson"])
          {
            string gjson_name = OUT_DIR+OUT_ID+"final_single_channel.geojson";
            source_points_data.print_data_to_geojson(gjson_name);
          }   
        }

        mod.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);
      }

      LSDRaster new_muddpile_topo_raster = mod.return_as_raster(); 
      if(this_bool_map["carve_before_fill"])
      {
        ct = new_muddpile_topo_raster.Breaching_Lindsay2016();
        ft = ct.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        ft = new_muddpile_topo_raster.fill(this_float_map["min_slope_for_fill"]);
      }

      // Write the raster
      string R_name = OUT_DIR+OUT_ID+"Final_DEM_fill";
      ft.write_raster(R_name,"bil");
      R_name = OUT_DIR+OUT_ID+"Final_DEM";
      new_muddpile_topo_raster.write_raster(R_name,"bil");
  
      cout << "\t Flow routing..." << endl;
      LSDFlowInfo FlowInfo(boundary_conditions,ft);

      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      cout << "\t Converting to flow area..." << endl;
      LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

      // calculate the distance from outlet
      cout << "\t Calculating flow distance..." << endl;
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

      cout << "\t Loading Sources..." << endl;
      cout << "\t Source file is... " << CHeads_file << endl;

      // load the sources
      vector<int> sources;
      if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
      {
        cout << endl << endl << endl << "==================================" << endl;
        cout << "The channel head file is null. " << endl;
        cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
        sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

        cout << "The number of sources is: " << sources.size() << endl;
      }
      else
      {
        cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
        sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
        cout << "\t Got sources!" << endl;
      }

      // now get the junction network
      LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

      // get the chi coordinate
      float A_0 = 1.0;
      float thresh_area_for_chi = float( this_int_map["threshold_contributing_pixels"] ); // needs to be smaller than threshold channel area

      float movern = this_float_map["m"]/this_float_map["n"];
      LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

      vector<int> BaseLevelJunctions = JunctionNetwork.get_BaseLevelJunctions();

      cout << endl << endl << endl << "=================" << endl;
      for (int i = 0; i< int(BaseLevelJunctions.size()); i++)
      {
        cout << "bl["<<i<<"]: " << BaseLevelJunctions[i] << endl;
      }
      cout <<"=================" << endl;

      vector<int> source_nodes;
      vector<int> outlet_nodes;
      vector<int> baselevel_node_of_each_basin;
      int n_nodes_to_visit = 10;
      JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

      LSDChiTools ChiTool_chi_checker(FlowInfo);
      ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                          ft, DistanceFromOutlet,
                          DrainageArea, chi_coordinate);

      string chi_data_maps_string = OUT_DIR+OUT_ID+"_final_chi_data_map.csv";
      ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_final_chi_data_map.geojson";
        LSDSpatialCSVReader thiscsv(chi_data_maps_string);
        thiscsv.print_data_to_geojson(gjson_name);
      }  
    }

    cout << "Thanks for snapping!" << endl;
  }

  //=======================================================================================
  //  /$$       /$$   /$$     /$$                                                        
  // | $$      |__/  | $$    | $$                                                        
  // | $$       /$$ /$$$$$$  | $$$$$$$   /$$$$$$   /$$$$$$$ /$$$$$$$   /$$$$$$   /$$$$$$ 
  // | $$      | $$|_  $$_/  | $$__  $$ /$$__  $$ /$$_____/| $$__  $$ |____  $$ /$$__  $$
  // | $$      | $$  | $$    | $$  \ $$| $$  \ $$|  $$$$$$ | $$  \ $$  /$$$$$$$| $$  \ $$
  // | $$      | $$  | $$ /$$| $$  | $$| $$  | $$ \____  $$| $$  | $$ /$$__  $$| $$  | $$
  // | $$$$$$$$| $$  |  $$$$/| $$  | $$|  $$$$$$/ /$$$$$$$/| $$  | $$|  $$$$$$$| $$$$$$$/
  // |________/|__/   \___/  |__/  |__/ \______/ |_______/ |__/  |__/ \_______/| $$____/ 
  //                                                                           | $$      
  //                                                                           | $$      
  //                                                                           |__/   
  // 
  //   /$$$$$$            /$$   /$$              /$$                              
  //  /$$__  $$          |__/  | $$             | $$                              
  // | $$  \__/  /$$$$$$  /$$ /$$$$$$   /$$$$$$$| $$  /$$$$$$   /$$$$$$   /$$$$$$ 
  // | $$       /$$__  $$| $$|_  $$_/  /$$_____/| $$ /$$__  $$ /$$__  $$ /$$__  $$
  // | $$      | $$  \__/| $$  | $$   |  $$$$$$ | $$| $$  \ $$| $$  \ $$| $$$$$$$$
  // | $$    $$| $$      | $$  | $$ /$$\____  $$| $$| $$  | $$| $$  | $$| $$_____/
  // |  $$$$$$/| $$      | $$  |  $$$$//$$$$$$$/| $$|  $$$$$$/| $$$$$$$/|  $$$$$$$
  //  \______/ |__/      |__/   \___/ |_______/ |__/ \______/ | $$____/  \_______/
  //                                                          | $$                
  //                                                          | $$                
  //                                                          |__/                
  //    
  //=======================================================================================
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
    // SMM: This is needed to get the initial lithological units for the first snap
    if (this_float_map["elevation_change"] != 0)
    {
      elev_raster.AdjustElevation(this_float_map["elevation_change"]);
      cout << "Adjusted the elevation of your raster by " << this_float_map["elevation_change"] << endl;
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
    if( this_bool_map["load_lithocode_to_K_csv"])
    {   
      string this_lookup_table_K = DATA_DIR+this_string_map["lookup_filename_K"];
      cout << "I am going to load a lookup table for K" << endl;
      cout << "The filename is: " << this_lookup_table_K << endl;
      LSDLC.read_strati_to_K_csv(this_lookup_table_K);
    }
    else
    {
      LSDLC.hardcode_fill_strati_to_K_map();  
    }  
    // Let's create the map of stratigraphy codes and K values    
    if( this_bool_map["load_lithocode_to_Sc_csv"])
    {   
      string this_lookup_table_Sc = DATA_DIR+this_string_map["lookup_filename_Sc"];
      cout << "I am going to load a lookup table for Sc" << endl;
      cout << "The filename is: " << this_lookup_table_Sc << endl;
      LSDLC.read_strati_to_Sc_csv(this_lookup_table_Sc);
    }
    else
    {
      LSDLC.hardcode_fill_strati_to_Sc_map();  
    }  
    

    // Now we'll get the lithology codes, convert them to K values and make a raster of K values for the model to use
    cout << "Getting the lithology codes" << endl;
    LSDIndexRaster lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(elev_raster,forbidden_lithocodes);
    cout << "Converting the litho codes to K values" << endl;
    LSDRaster K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);
    cout << "Converting the litho codes to Sc values" << endl;
    LSDRaster Sc_values = LSDLC.index_to_Sc(lithocodes_index_raster, 0.21);
        
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
    // LSDSpatialCSVReader single_channel_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

    if (this_float_map["single_channel_drop"] != 0)
    {
      cout << endl << endl << endl << "==================================" << endl;
      cout << "I am dropping and adjusting the slope of the single channel." << endl;
      cout << "The data keys are: " << endl;
      source_points_data.print_data_map_keys_to_screen();
      cout << "==================================" << endl << endl << endl <<endl;

      string elevation_column_name = this_string_map["single_channel_elev_string"];
      string fd_column_name = this_string_map["single_channel_fd_string"];

      source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
      source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]);
    }


    if (this_bool_map["only_test_impose_single_channel"])
    {
      cout << "=====================================" << endl;
      cout << "I'm not going to snap, I'm just imposing the channels" << endl;
      string csv_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.csv";
      source_points_data.print_data_to_csv(csv_name);
      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"testing_single_channel_imposed.geojson";
        source_points_data.print_data_to_geojson(gjson_name);
      }   

      LSDRaster testing_topo_raster = mod.return_as_raster();
      LSDRasterMaker RM(testing_topo_raster);
      RM.impose_channels_with_buffer(source_points_data, this_float_map["min_slope_for_fill"]*2, this_string_map["single_channel_elev_string"]);

      LSDRaster imposed_raster = RM.return_as_raster();
      string irastername = OUT_DIR+OUT_ID+"_imposed";
      imposed_raster.write_raster(irastername,"bil");

      exit(0);
    }


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
      cout << "The U value at (500,500) is: " << U_values.get_data_element(500,500) << endl;

      cout << "*****************************************" << endl;
      cout << "This is snap " << i+1 << " of " << n_snaps << endl;
      cout << "*****************************************" << endl;
      


      // if(i>0)
      // {
      cout << "Let me create a raster of the topography from the previous loop." << endl;
      new_muddpile_topo_raster = mod.return_as_raster();
      cout << "Now I'll check where in the lithocube we are and make a raster of lithocodes." << endl;
      lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(new_muddpile_topo_raster,forbidden_lithocodes);
      cout << "I am now converting the lithocodes index raster to an LSDRaster" << endl;
      lithocodes_raster = LSDRaster(lithocodes_index_raster);
      cout << "Converting the litho codes to K values" << endl;
      K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);
      cout << "Converting the litho codes to Sc values" << endl;
      Sc_values = LSDLC.index_to_Sc(lithocodes_index_raster, 0.21);

      // }

      mod.fluvial_snap_to_steady_variable_K_variable_U(K_values, U_values, source_points_data,this_bool_map["carve_before_fill"], this_string_map["single_channel_elev_string"]);


      if (this_bool_map["cycle_hillslope_snapping"])
      {
        mod.hillslope_snap_to_steady_variable_Sc(Sc_values, this_bool_map["carve_before_fill"], this_int_map["threshold_contributing_pixels"]);
      }

      if (this_bool_map["use_initial_topography_as_max_elevation"])
      {
        mod.cap_elevations(InitialRaster);
      }

      // Snap to critical slopes now? 
      if (i == n_snaps-1)
      {
        cout << "I've got to the last snap and need to do some special stuff" << endl;
        // first, get the channels for burning
        // make sure the single channel is imposed
        mod.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);

        LSDSpatialCSVReader chan_data = mod.get_channels_for_burning(this_int_map["threshold_contributing_pixels"],
                                                                     path_name,this_string_map["temp_chan_name_for_crit_slopes"]);


        if (this_bool_map["cycle_hillslope_snapping"])
        {
          cout << "I snapped your critical slopes during each snapping cycle" << endl;
          critical_slopes = mod.return_as_raster();
        }
        else
        {
          cout << "Last snap! Getting critical slopes and smoothing" << endl;
          critical_slopes = mod.basic_valley_fill_critical_slope(Sc_values,this_int_map["threshold_contributing_pixels"]);
        }
        


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
        RM.impose_channels(chan_data, this_string_map["single_channel_elev_string"]);
        RM.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);
        

        if (this_bool_map["use_initial_topography_as_max_elevation"])
        {
          mod.cap_elevations(InitialRaster);
          RM.cap_elevations(InitialRaster);
        }
        critical_slopes = RM.return_as_raster();

        cout << "Final critical slope snapping, I will give you the lithocodes, K and Sc." << endl;
        cout << "Now I'll check where in the lithocube we are and make a raster of lithocodes." << endl;
        lithocodes_index_raster = LSDLC.get_lithology_codes_from_raster(critical_slopes,forbidden_lithocodes);
        cout << "I am now converting the lithocodes index raster to an LSDRaster" << endl;
        lithocodes_raster = LSDRaster(lithocodes_index_raster);
        cout << "Converting the litho codes to K values" << endl;
        K_values = LSDLC.index_to_K(lithocodes_index_raster, 0.000003);
        cout << "Converting the litho codes to Sc values" << endl;
        Sc_values = LSDLC.index_to_Sc(lithocodes_index_raster, 0.21);
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


    if (this_bool_map["print_final_channel_network"])
    {
      LSDRaster ft,ct;

      if(this_bool_map["impose_single_channel_printing"])
      {
        // load the single channel
        LSDRasterInfo RI(InitialRaster);
        // Get the latitude and longitude
        //cout << "I am reading points from the file: "+ this_string_map["channel_source_fname"] << endl;
        LSDSpatialCSVReader source_points_data( RI, (DATA_DIR+this_string_map["fixed_channel_csv_name"]) );

        if(this_bool_map["fixed_channel_dig_and_slope"])
        {
          string fd_column_name = this_string_map["single_channel_fd_string"];
          string elevation_column_name = this_string_map["single_channel_elev_string"];

          source_points_data.data_column_add_float(elevation_column_name, -this_float_map["single_channel_drop"]);
          source_points_data.enforce_slope(fd_column_name, elevation_column_name, this_float_map["min_slope_for_fill"]); 

          string csv_name = OUT_DIR+OUT_ID+"final_single_channel.csv";
          source_points_data.print_data_to_csv(csv_name);
          if ( this_bool_map["convert_csv_to_geojson"])
          {
            string gjson_name = OUT_DIR+OUT_ID+"final_single_channel.geojson";
            source_points_data.print_data_to_geojson(gjson_name);
          }   
        }

        mod.impose_channels(source_points_data, this_string_map["single_channel_elev_string"]);
      }

      LSDRaster new_muddpile_topo_raster = mod.return_as_raster(); 
      if(this_bool_map["carve_before_fill"])
      {
        ct = new_muddpile_topo_raster.Breaching_Lindsay2016();
        ft = ct.fill(this_float_map["min_slope_for_fill"]);
      }
      else
      {
        ft = new_muddpile_topo_raster.fill(this_float_map["min_slope_for_fill"]);
      }

      // Write the raster
      string R_name = OUT_DIR+OUT_ID+"Final_DEM_fill";
      ft.write_raster(R_name,"bil");
      R_name = OUT_DIR+OUT_ID+"Final_DEM";
      new_muddpile_topo_raster.write_raster(R_name,"bil");

      cout << "\t Flow routing..." << endl;
      LSDFlowInfo FlowInfo(boundary_conditions,ft);

      // calculate the flow accumulation
      cout << "\t Calculating flow accumulation (in pixels)..." << endl;
      LSDIndexRaster FlowAcc = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();

      cout << "\t Converting to flow area..." << endl;
      LSDRaster DrainageArea = FlowInfo.write_DrainageArea_to_LSDRaster();

      // calculate the distance from outlet
      cout << "\t Calculating flow distance..." << endl;
      LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

      cout << "\t Loading Sources..." << endl;
      cout << "\t Source file is... " << CHeads_file << endl;

      // load the sources
      vector<int> sources;
      if (CHeads_file == "NULL" || CHeads_file == "Null" || CHeads_file == "null")
      {
        cout << endl << endl << endl << "==================================" << endl;
        cout << "The channel head file is null. " << endl;
        cout << "Getting sources from a threshold of "<< this_int_map["threshold_contributing_pixels"] << " pixels." <<endl;
        sources = FlowInfo.get_sources_index_threshold(FlowAcc, this_int_map["threshold_contributing_pixels"]);

        cout << "The number of sources is: " << sources.size() << endl;
      }
      else
      {
        cout << "Loading channel heads from the file: " << DATA_DIR+CHeads_file << endl;
        sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv",2);
        cout << "\t Got sources!" << endl;
      }

      // now get the junction network
      LSDJunctionNetwork JunctionNetwork(sources, FlowInfo);

      // get the chi coordinate
      float A_0 = 1.0;
      float thresh_area_for_chi = float( this_int_map["threshold_contributing_pixels"] ); // needs to be smaller than threshold channel area

      float movern = this_float_map["m"]/this_float_map["n"];
      LSDRaster chi_coordinate = FlowInfo.get_upslope_chi_from_all_baselevel_nodes(movern,A_0,thresh_area_for_chi);

      vector<int> BaseLevelJunctions = JunctionNetwork.get_BaseLevelJunctions();

      cout << endl << endl << endl << "=================" << endl;
      for (int i = 0; i< int(BaseLevelJunctions.size()); i++)
      {
        cout << "bl["<<i<<"]: " << BaseLevelJunctions[i] << endl;
      }
      cout <<"=================" << endl;

      vector<int> source_nodes;
      vector<int> outlet_nodes;
      vector<int> baselevel_node_of_each_basin;
      int n_nodes_to_visit = 10;
      JunctionNetwork.get_overlapping_channels_to_downstream_outlets(FlowInfo, BaseLevelJunctions, DistanceFromOutlet,
                                  source_nodes,outlet_nodes,baselevel_node_of_each_basin,n_nodes_to_visit);

      LSDChiTools ChiTool_chi_checker(FlowInfo);
      ChiTool_chi_checker.chi_map_automator_chi_only(FlowInfo, source_nodes, outlet_nodes, baselevel_node_of_each_basin,
                          ft, DistanceFromOutlet,
                          DrainageArea, chi_coordinate);

      string chi_data_maps_string = OUT_DIR+OUT_ID+"_final_chi_data_map.csv";
      ChiTool_chi_checker.print_chi_data_map_to_csv(FlowInfo, chi_data_maps_string);

      if ( this_bool_map["convert_csv_to_geojson"])
      {
        string gjson_name = OUT_DIR+OUT_ID+"_final_chi_data_map.geojson";
        LSDSpatialCSVReader thiscsv(chi_data_maps_string);
        thiscsv.print_data_to_geojson(gjson_name);
      }  
    }
  }


  //==============================================================================================
  //
  //  /$$$$$$$$                                     /$$                                        
  // |__  $$__/                                    |__/                                        
  //    | $$  /$$$$$$  /$$$$$$  /$$$$$$$   /$$$$$$$ /$$  /$$$$$$  /$$$$$$$   /$$$$$$$  /$$$$$$ 
  //    | $$ /$$__  $$|____  $$| $$__  $$ /$$_____/| $$ /$$__  $$| $$__  $$ /$$_____/ /$$__  $$
  //    | $$| $$  \__/ /$$$$$$$| $$  \ $$|  $$$$$$ | $$| $$$$$$$$| $$  \ $$| $$      | $$$$$$$$
  //    | $$| $$      /$$__  $$| $$  | $$ \____  $$| $$| $$_____/| $$  | $$| $$      | $$_____/
  //    | $$| $$     |  $$$$$$$| $$  | $$ /$$$$$$$/| $$|  $$$$$$$| $$  | $$|  $$$$$$$|  $$$$$$$
  //    |__/|__/      \_______/|__/  |__/|_______/ |__/ \_______/|__/  |__/ \_______/ \_______/
  // 
  //============================================================================================     
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
    //mod.set_float_print_interval(this_float_map["float_print_interval"]);

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

  if( this_bool_map["test_transient_channel"])
  {
    // get the single channel
    cout << "Let me get the single channel...";
    LSDRasterInfo RI(InitialRaster);
    LSDSpatialCSVReader transient_channel_data( RI, (DATA_DIR+this_string_map["transient_channel_csv_name"]) );
    cout << "Got it!" << endl;

    string fd_column_name = this_string_map["single_channel_fd_string"];
    string elevation_column_name = this_string_map["single_channel_elev_string"];


    vector<float> phases;
    vector< vector<float> > elevations_vecvec;
    mod.process_baselevel( transient_channel_data, this_string_map["transient_channel_timing_prefix"], 
                       this_float_map["transient_channel_timing_multiplier"], this_int_map["transient_channel_n_time_columns"], 
                       this_int_map["transient_channel_phase_steps"],
                       phases,elevations_vecvec);

    exit(0);
  }


  //========================================================================================================================
  //
  //.######..#####....####...##..##...####...######..######..##..##..######......##......######..######..##..##...####..
  //...##....##..##..##..##..###.##..##........##....##......###.##....##........##........##......##....##..##..##..##.
  //...##....#####...######..##.###...####.....##....####....##.###....##........##........##......##....######..##..##.
  //...##....##..##..##..##..##..##......##....##....##......##..##....##........##........##......##....##..##..##..##.
  //...##....##..##..##..##..##..##...####...######..######..##..##....##........######..######....##....##..##...####..
  //
  //========================================================================================================================  
  if( this_bool_map["run_transient_base_level"] )
  {
    cout << "Hi there. Let me test some transient runs for you. " << endl;

    cout << endl << endl << endl << endl << "Before I start, I need to make sure the raster is filled. " << endl;
    LSDRaster temp_raster = mod.return_as_raster();
    LSDRaster carved_topography;
    LSDRaster filled_topography;

    if(this_bool_map["carve_before_fill"])
    {
      carved_topography = temp_raster.Breaching_Lindsay2016();
      filled_topography = carved_topography.fill(this_float_map["min_slope_for_fill"]);
    }
    else
    {
      filled_topography = temp_raster.fill(this_float_map["min_slope_for_fill"]);
    }
    LSDRasterModel temp_mod2(filled_topography);
    mod = temp_mod2;
    cout << "Done filling raster" << endl;

    cout << "Let me check the baselevel situation" << endl;
    LSDFlowInfo FI(filled_topography);
    cout << "N baselevel nodes is " << FI.get_NBaseLevelNodes() << endl << endl << endl;



    // set the print interval
    mod.set_print_interval( this_int_map["print_interval"] );

    // set the soil transport coefficient
    mod.set_D( this_float_map["D"]);

    cout << "The printing interval is: " << mod.get_float_print_interval() << " years." << endl;

    // set the end time
    current_end_time = this_float_map["transient_maximum_time"];
    mod.set_endTime(current_end_time);

    // Set the timestep
    cout << "The timestep is " << this_float_map["dt"] << endl;
    mod.set_timeStep( this_float_map["dt"] );    


    LSDLithoCube LSDLC; 
    LSDRasterInfo RI(InitialRaster);
  
    cout <<"Let me get some K, Sc and U rasters for you." << endl;
    cout << "If you are using the lithocube the K and Sc rasters will be ignored." << endl;
    string U_name = OUT_DIR+OUT_ID+"_ConstURaster";
    string rext = "bil";
    LSDRaster KRaster;
    LSDRaster ScRaster;

    cout << "Getting the constant uplift raster...";
    LSDRaster URaster(U_name,rext);;
    cout << "got it" << endl;


    // First check to see if you are using the lithocube
    if (this_bool_map["use_lithocube_for_transient"])
    {
      cout << "  I am going to use a lithocube in your transient runs. " << endl;
      
      // Then load the vo and asc files containing the 3D litho data 
      string this_vo_filename = DATA_DIR+this_string_map["vo_filename"];
      string this_vo_asc = DATA_DIR+this_string_map["vo_asc_filename"];

      // Set names for K raster and lithocode raster to be printed 
      string K_name = DATA_DIR+OUT_ID+"_KRasterLSDLC";
      string Sc_name = DATA_DIR+OUT_ID+"_ScRasterLSDLC";
      string lithocodes_raster_name = DATA_DIR+OUT_ID+"_LithoCodes";

      
      // Initialise the lithocube and ingest the lithology data
      cout << "I am going to load lithology information" << endl;
      cout << "The filename is: " << this_vo_filename << endl;
      LSDLithoCube LSDLC_ingest(this_vo_filename);
      LSDLC_ingest.ingest_litho_data(this_vo_asc);   

      // Deal with the georeferencing
      int UTM_zone;
      bool is_North = true;
      string NorS = "N";
      RI.get_UTM_information(UTM_zone,is_North);
      if (is_North == false)
      {
        NorS = "S";
      }  
      LSDLC_ingest.impose_georeferencing_UTM(UTM_zone, NorS);

      // Let's create the maps of stratigraphy codes and K values, and Sc values 
      if( this_bool_map["load_lithocode_to_K_csv"])
      {   
        string this_lookup_table_K = DATA_DIR+this_string_map["lookup_filename_K"];
        cout << "I am going to load a lookup table for K" << endl;
        cout << "The filename is: " << this_lookup_table_K << endl;
        LSDLC_ingest.read_strati_to_K_csv(this_lookup_table_K);
      }
      else
      {
        LSDLC_ingest.hardcode_fill_strati_to_K_map();  
      }  
      // Let's create the map of stratigraphy codes and K values    
      if( this_bool_map["load_lithocode_to_Sc_csv"])
      {   
        string this_lookup_table_Sc = DATA_DIR+this_string_map["lookup_filename_Sc"];
        cout << "I am going to load a lookup table for Sc" << endl;
        cout << "The filename is: " << this_lookup_table_Sc << endl;
        LSDLC_ingest.read_strati_to_Sc_csv(this_lookup_table_Sc);
      }
      else
      {
        LSDLC_ingest.hardcode_fill_strati_to_Sc_map();  
      }  
      
      // Now we'll get the lithology codes, convert them to K values and make a raster of K values for the model to use
      cout << "Getting the lithology codes" << endl;
      LSDIndexRaster lithocodes_index_raster = LSDLC_ingest.get_lithology_codes_from_raster(InitialRaster,forbidden_lithocodes);
      cout << "Converting the litho codes to K values" << endl;
      LSDRaster K_values = LSDLC_ingest.index_to_K(lithocodes_index_raster, 0.000003);
      cout << "Converting the litho codes to Sc values" << endl;
      LSDRaster Sc_values = LSDLC_ingest.index_to_Sc(lithocodes_index_raster, 0.21);

      KRaster = K_values;
      ScRaster = Sc_values;     
      LSDLC = LSDLC_ingest;
    }
    else
    {
      cout << "  This transient run will use a static K raster (not a lithocube that will change K values as it is exhumed)." << endl;
      cout << "  Let me read the rasters for you" << endl;
      string K_name = OUT_DIR+OUT_ID+"_ConstKRaster";
      string Sc_name = OUT_DIR+OUT_ID+"_ConstScRaster";
      LSDRaster K_values(K_name,rext);
      LSDRaster Sc_values(Sc_name,rext);
      KRaster = K_values;
      ScRaster = Sc_values;  

      cout << "  Got the constant rasters. " << endl;
    }


    // get the single channel
    cout << "Let me get the single channel...";
    LSDSpatialCSVReader transient_channel_data( RI, (DATA_DIR+this_string_map["transient_channel_csv_name"]) );
    cout << "Got it!" << endl;

    string fd_column_name = this_string_map["single_channel_fd_string"];
    string elevation_column_name = this_string_map["single_channel_elev_string"];

    // This step takes the data that has elevations through timesteps and 
    // extracts the elevations as a function of times (which we call "phases")
    cout << "Now I need to process the baselevel information." << endl;
    cout << "WARNING: If I have a baselevel code for reading area, the area column needs to be: area"  << endl;
    cout << "At the moment you can only change this in the source code (LSDRasterModel line 7978 or thereabouts" << endl;
    vector<float> phases;
    vector< vector<float> > elevations_vecvec;
    vector<float> outlet_elevations;

    if (this_bool_map["use_transient_outlet"])
    {
      cout << "I am processing a transient run with only a single outlet." << endl;
      cout << "You need a transient csv for this. " << endl;
      string transient_infile_name = DATA_DIR+this_string_map["transient_outlet_fname"];
      cout << "The name of the transient infile is: " << transient_infile_name << endl;
      mod.process_transient_file(transient_infile_name, phases,outlet_elevations);
    }
    else
    {
      cout << "I am processing baselevel with a full transient channel" << endl;
      cout << "This has different elevations at pixels throughout" << endl;
      cout << "one or more channels." << endl;
      mod.process_baselevel( transient_channel_data, this_string_map["transient_channel_timing_prefix"], 
                       this_float_map["transient_channel_timing_multiplier"], this_int_map["transient_channel_n_time_columns"], 
                       this_int_map["transient_channel_phase_steps"],
                       phases,elevations_vecvec);
    }
    cout << "Finished with the baselevel" << endl;



    // Make sure the baselevel channel does not sit above the rest of the raster
    if (this_bool_map["raise_raster_over_base_level"])
    {
      cout << "Min an max elevations are: " << mod.min_elevation() << " " << mod.max_elevation() << endl;
      float min_elevation = mod.min_elevation();
      
      vector<float> this_elev = transient_channel_data.data_column_to_float(this_string_map["single_channel_elev_string"]);

      float outlet = this_elev[0];
      cout << "Outlet elevation is: " << outlet << endl;
      if (outlet > min_elevation)
      {
        mod.add_fixed_elevation( outlet-min_elevation+2);
      }
    }

    // Add the baselevel switch to the baselevel channel
    mod.baselevel_run_area_switch( transient_channel_data, this_string_map["baselevel_switch_column"]);

    cout << "I am going to turn off nonlinear diffusion." << endl;
    cout << "In this model the hillslopes are critical slopes near channel." << endl;
    cout << "And linear diffusion on the ridgetops." << endl;
    mod.set_hillslope(false);

    // Make a cumulative uplift raster  
    // We get the right shape and then populate with a single value of 0 
    cout << "Making the cumulative uplift raster." << endl;
    LSDRaster Cumulative_uplift = mod.return_as_raster();
    Cumulative_uplift = Cumulative_uplift.PopulateRasterSingleValue(0.0);
    Cumulative_uplift.write_raster(OUT_DIR+OUT_ID+"_test_CU_","bil");
    
    cout << "Let me print the rasters before I do any simulation." << endl;
    int this_frame = 9999;
    mod.print_rasters( this_frame );    

    // reset the frame
    mod.set_current_frame(0);

    // Set the timestep
    cout << "I am repeating the timestep enforcement: " << this_float_map["dt"] << endl;
    mod.set_timeStep( this_float_map["dt"] );    

    cout << "Okay, now into the simulation loop...";
    if (this_bool_map["use_adaptive_timestep"])
    {
      cout << "Using an adaptive timestep. This might slow down your model a lot (but it ensures it is stable)! "<< endl;
    }

    if (this_bool_map["use_lithocube_for_transient"]) 
    {
      cout << endl << endl << "=============================================" << endl;
      cout << "I repeat, I am running the transient model using the lithocube" << endl;

      if (this_bool_map["use_transient_outlet"])
      {
        cout << "I am running the model with a lowering outlet node" << endl;
        mod.run_components_combined_imposed_baselevel( URaster, this_bool_map["use_adaptive_timestep"], transient_channel_data,
                                                  phases,outlet_elevations, Cumulative_uplift, LSDLC, forbidden_lithocodes,
                                                  this_bool_map["print_lithocode_raster"], this_bool_map["use_hillslope_hybrid"],
                                                  this_int_map["threshold_contributing_pixels"], 
                                                  this_float_map["min_slope_for_fill"],
                                                  this_bool_map["print_exhumation_and_cumulative_uplift"], this_string_map["single_channel_elev_string"]);        
      }
      else
      {
        cout << "I am running the model with an imposed channel" << endl;
        mod.run_components_combined_imposed_baselevel( URaster, this_bool_map["use_adaptive_timestep"], transient_channel_data,
                                                  phases,elevations_vecvec, Cumulative_uplift, LSDLC, forbidden_lithocodes,
                                                  this_bool_map["print_lithocode_raster"], this_bool_map["use_hillslope_hybrid"],
                                                  this_int_map["threshold_contributing_pixels"], 
                                                  this_float_map["min_slope_for_fill"],
                                                  this_bool_map["print_exhumation_and_cumulative_uplift"], this_string_map["single_channel_elev_string"]);
      }
           
    } 
    else
    {
      cout << endl << endl << "=============================================" << endl;
      cout << "I repeat, I am running the transient model with a fixed K raster" << endl;
      mod.run_components_combined_imposed_baselevel( URaster, KRaster, this_bool_map["use_adaptive_timestep"], transient_channel_data,
                                                  phases,elevations_vecvec, this_float_map["min_slope_for_fill"], this_string_map["single_channel_elev_string"]);

    }
      
    cout << "Done!" << endl;
 
  }

}
