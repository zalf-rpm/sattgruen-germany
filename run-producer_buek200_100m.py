#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import time
import os
import math
import json
import csv
import copy
#from StringIO import StringIO
from datetime import date, timedelta
from collections import defaultdict
import sys
import zmq

import sqlite3
import sqlite3 as cas_sq3
import numpy as np
from pyproj import Proj, transform, CRS, Transformer
import PIL, PIL.Image as Image
PIL.Image.MAX_IMAGE_PIXELS = 6940000000
#from PIL import TiffTags

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib


# chose setup for running as container (in docker) or local run 
# for local run you need a local monica server running e.g "monica-zmq-server -bi -i tcp://*:6666 -bo -o tcp://*:7777"
# if you use a local monica-zmq-server, set "monica-path-to-climate-dir" to the folder where your climate data is found

# or a local docker image:
# docker run -p 6666:6666 -p 7777:7777 --env monica_instances=3 --rm --name monica_test -v climate-data:/monica_data/climate-data zalfrpm/monica-cluster:latest
# if you use docker, set "monica-path-to-climate-dir" = "/monica_data/climate-data/" 
# and create a volume for the climate data, e.g for a network drive
# docker volume create --driver local \
#     --opt type=cifs \fwat
#     --opt device='//network_drive_ip/archiv-daten/md/data/climate/' \
#     --opt o='username=your_username,password=your_password' \
# climate-data

PATHS = {
    # adjust the local path to your environment
    "cs-local-remote": {
        "include-file-base-path": "D:/awork/zalf/monica/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "D:/awork/zalf/monica/monica_example/monica-data/climate-data/", # mounted path to archive or hard drive with climate data 
        "monica-path-to-climate-dir": "/monica_data/climate-data/", # mounted path to archive accessable by monica executable
        "path-to-data-dir": "D:/awork/zalf/monica/monica_example/monica-data/data/", # mounted path to archive or hard drive with data 
        "path-to-projects-dir": "D:/awork/zalf/monica/monica_example/monica-data/projects/", # mounted path to archive or hard drive with project data 
    },
    "mbm-local-remote": {
        "include-file-base-path": "C:/Users/berg.ZALF-AD/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "W:/FOR/FPM/data/climate/", # mounted path to archive or hard drive with climate data 
        "monica-path-to-climate-dir": "/monica_data/climate-data/", # mounted path to archive accessable by monica executable
        "path-to-data-dir": "C:/Users/berg.ZALF-AD/GitHub/sattgruen-germany/monica-data/data/" # mounted path to archive or hard drive with data 
    },
    "hpc-remote": {
        "include-file-base-path": "/beegfs/common/GitHub/zalf-rpm/monica-parameters/",
        "path-to-climate-dir": "/beegfs/common/data/climate/", 
        "monica-path-to-climate-dir": "/monica_data/climate-data/", 
        "path-to-data-dir": "/beegfs/common/data/" 
    },
    "container": {
        "include-file-base-path": "/home/monica-parameters/", # monica parameter location in docker image
        "monica-path-to-climate-dir": "/monica_data/climate-data/",  # mounted path to archive on cluster docker image 
        "path-to-climate-dir": "/monica_data/climate-data/", # needs to be mounted there
        "path-to-data-dir": "/monica_data/data/", # needs to be mounted there
        "path-to-projects-dir": "/monica_data/project/", # needs to be mounted there
    },
    "remoteProducer-remoteMonica": {
        "include-file-base-path": "/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/data/", # mounted path to archive or hard drive with climate data 
        "monica-path-to-climate-dir": "/monica_data/climate-data/", # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./monica-data/data/" # mounted path to archive or hard drive with data 
    }
}

DEFAULT_HOST = "localhost"
DEFAULT_PORT = "6669"
RUN_SETUP = "[1]"
SETUP_FILE = "sim_setups_cut_grids.csv"
DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_HEIGHT = "germany/dem_1000_gk5.asc" 
DATA_GRID_SLOPE = "germany/slope_1000_gk5.asc"
DATA_GRID_LAND_USE = "germany/corine2006_1000_gk5.asc"
DATA_GRID_SOIL = "germany/BUEK200_1000_gk5.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}{climate_data}/csvs/latlon-to-rowcol.json"
TEMPLATE_PATH_CLIMATE_CSV = "{climate_data}/csvs/{climate_model_folder}{climate_scenario_folder}{climate_region}/row-{crow}/col-{ccol}.csv"
GEO_TARGET_GRID="epsg:31469" #proj4 -> 3-degree gauss-kruger zone 5 (=Germany) https://epsg.io/31469

DEBUG_DONOT_SEND = False
DEBUG_WRITE = False
DEBUG_ROWS = 10
DEBUG_WRITE_FOLDER = "./debug_out"
DEBUG_WRITE_CLIMATE = True

# some values in these templates will be overwritten by the setup 
TEMPLATE_SIM_JSON="sim.json" 
TEMPLATE_CROP_JSON="crop_template.json"
TEMPLATE_SITE_JSON="site.json"

# commandline parameters e.g "server=localhost port=6666 shared_id=2"
def run_producer(server = {"server": None, "port": None}, shared_id = None):
    "main"

    context = zmq.Context()
    socket = context.socket(zmq.PUSH)
    #config_and_no_data_socket = context.socket(zmq.PUSH)

    config = {
        "mode": "remoteProducer-remoteMonica",
        "server-port": server["port"] if server["port"] else DEFAULT_PORT,
        "server": server["server"] if server["server"] else DEFAULT_HOST,
        "start-row": "0", 
        "end-row": "-1",
        "sim.json": TEMPLATE_SIM_JSON,
        "crop.json": TEMPLATE_CROP_JSON,
        "site.json": TEMPLATE_SITE_JSON,
        "setups-file": SETUP_FILE,
        "run-setups": RUN_SETUP,
        "shared_id": shared_id
    }
    
    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print("config:", config)

    # select paths 
    paths = PATHS[config["mode"]]
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)
    #soil_db_con = cas_sq3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB) #CAS.
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"], flush=True)

    #transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    wgs84 = Proj(init="epsg:4326") #proj4 -> (World Geodetic System 1984 https://epsg.io/4326)
    gk5 = Proj(init=GEO_TARGET_GRID) 
    #etrs89_crs = CRS.from_epsg(3035)
    latlon_crs = CRS.from_epsg(4326)
    gk3_crs = CRS.from_epsg(31467)
    gk5_crs = CRS.from_epsg(31469)
    #latlon_to_gk3 = Transformer.from_crs(latlon_crs, gk3_crs, always_xy=True)
    gk5_to_gk3 = Transformer.from_crs(gk5_crs, gk3_crs, always_xy=True)
    gk3_to_gk5 = Transformer.from_crs(gk3_crs, gk5_crs, always_xy=True)
    #latlon_to_etrs89 = Transformer.from_crs(latlon_crs, etrs89_crs, always_xy=True)
    #im = Image.open("Tif_to_Asc_Mowing/mowingEvents_PYP_2020_sum_sieve_4_50_na.tif")
    #print(im.getpixel((30000, 30000)))

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2

    # cutting grids
    cutting_year_start = 2017
    cutting_year_end = 2019
    band_start = 0
    band_end = 9
    mow_grid_subfilepath = {
            "no_of_cuts": "{year}/mowingEvents_PYP_{year}_sum_sieve_4_50_na_31467_rsmpl_100.asc.gz",
            "doys": "{year}/mowing_DOY_PYP_{year}DOY_{band}_31467_rsmpl_100.asc.gz"
    }
    mow_grids = {
        "no_of_cuts": {},
        "doys": defaultdict(list)
    }
    path_to_mow_grid_dir = paths["path-to-climate-dir"] + "management/sattgruen_cutting_dates/"
    mow_metadata, _ = Mrunlib.read_header(path_to_mow_grid_dir + mow_grid_subfilepath["no_of_cuts"].format(year=cutting_year_start))
    for year in range(cutting_year_start, cutting_year_end+1):
        mow_grids["no_of_cuts"][year] = np.loadtxt(path_to_mow_grid_dir + mow_grid_subfilepath["no_of_cuts"].format(year=year), dtype=int, skiprows=6)
        print("no_cuts mow_grid loaded for", year, flush=True)
        for band in range(band_start, band_end+1):
            mow_grids["doys"][year].append(np.loadtxt(path_to_mow_grid_dir + mow_grid_subfilepath["doys"].format(year=year, band=band), dtype=int, skiprows=6))
            print("doys mow_grid loaded for", year, "and band", band, flush=True)
    #mow_gk3_interpolate = Mrunlib.create_ascii_grid_interpolator(mow_grids["no_of_cuts"][cutting_year_start], mow_metadata, return_row_col=True)
    print("mow_grid row col interpolator created", flush=True)
    
    # height data for germany
    path_to_dem_grid = paths["path-to-data-dir"] + DATA_GRID_HEIGHT 
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=int, skiprows=6)
    dem_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid, flush=True)
    
    # slope data
    path_to_slope_grid = paths["path-to-data-dir"] + DATA_GRID_SLOPE
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid, flush=True)

    # land use data
    path_to_corine_grid = paths["path-to-data-dir"] + DATA_GRID_LAND_USE
    corine_meta, _ = Mrunlib.read_header(path_to_corine_grid)
    corine_grid = np.loadtxt(path_to_corine_grid, dtype=int, skiprows=6)
    corine_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(corine_grid, corine_meta)
    print("read: ", path_to_corine_grid, flush=True)

    # soil data
    path_to_soil_grid = paths["path-to-data-dir"] + DATA_GRID_SOIL
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    soil_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(soil_grid, soil_metadata)
    print("read: ", path_to_soil_grid, flush=True)

    cdict = {}
    climate_data_to_gk5_interpolator = {}
    for run_id in run_setups:
        setup = setups[run_id]
        climate_data = setup["climate_data"]
        if not climate_data in climate_data_to_gk5_interpolator:
            # path to latlon-to-rowcol.json
            path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"], climate_data=climate_data)
            climate_data_to_gk5_interpolator[climate_data] = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84, gk5, cdict)
            print("created climate_data to gk5 interpolator: ", path, flush=True)

    sent_env_count = 1
    start_time = time.clock()

    listOfClimateFiles = set()
    # run calculations for each setup
    for _, setup_id in enumerate(run_setups):

        if setup_id not in setups:
            continue
        start_setup_time = time.clock()      

        setup = setups[setup_id]
        climate_data = setup["climate_data"]
        climate_model = setup["climate_model"]
        climate_scenario = setup["climate_scenario"]
        climate_region = setup["climate_region"]

        # read template sim.json 
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date acording to setup
        if setup["start_year"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_year"]) + "-01-01"
        if setup["end_year"]:
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_year"]) + "-12-31" 
        sim_json["include-file-base-path"] = paths["include-file-base-path"]

        # read template site.json 
        with open(setup.get("site.json", config["site.json"])) as _:
            site_json = json.load(_)
        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)
        crop_json["cropRotation"][0]["worksteps"][0]["date"] = crop_json["cropRotation"][0]["worksteps"][0]["date"].replace("XXXX", str(setup["start_year"]))

        # create environment template from json templates
        env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })
        # set shared id in template
        if config["shared_id"]:
            env_template["sharedId"] = config["shared_id"]

        base_cols = int(mow_metadata["ncols"])
        base_rows = int(mow_metadata["nrows"])
        scellsize = int(mow_metadata["cellsize"])
        xllcorner = int(mow_metadata["xllcorner"])
        yllcorner = int(mow_metadata["yllcorner"])

        cutting_cm_template = env_template["cropRotation"].pop(1)
        cutting_ws_template = cutting_cm_template["worksteps"].pop(1)
        crop_rotation_template = env_template.pop("cropRotation")

        #unknown_soil_ids = set()
        soil_id_cache = {}
        #print("All Rows x Cols: " + str(base_rows) + "x" + str(base_cols), flush=True)
        base_grid = mow_grids["no_of_cuts"][2017]
        for brow in range(0, base_rows):
            #print(brow,)

            if brow < int(config["start-row"]):
                continue
            elif int(config["end-row"]) > 0 and brow > int(config["end-row"]):
                break

            #if brow != 102:
            #    continue

            for bcol in range(0, base_cols):

            #    if bcol != 258:
            #        continue

                if base_grid[brow, bcol] == -9999:
                    continue

                bh_gk3 = yllcorner + (scellsize / 2) + (base_rows - brow - 1) * scellsize
                br_gk3 = xllcorner + (scellsize / 2) + bcol * scellsize
                br_gk5, bh_gk5 = gk3_to_gk5.transform(br_gk3, bh_gk3)
                         
                soil_id = int(soil_gk5_interpolate(br_gk5, bh_gk5))
                if soil_id in soil_id_cache:
                    sp_json = soil_id_cache[soil_id]
                else:
                    sp_json = soil_io3.soil_parameters(soil_db_con, soil_id)
                    soil_id_cache[soil_id] = sp_json

                if len(sp_json) == 0:
                    print("row/col:", brow, "/", bcol, "has unknown soil_id:", soil_id)
                    #unknown_soil_ids.add(soil_id)
                    continue
               
                #inter = crow/ccol encoded into integer
                crow, ccol = climate_data_to_gk5_interpolator[climate_data](br_gk5, bh_gk5)

                # check if current grid cell is used for agriculture                
                if setup["landcover"]:
                    corine_id = corine_gk5_interpolate(br_gk5, bh_gk5)
                    if corine_id not in [200, 210, 211, 212, 240, 241, 242, 243, 244]:
                        continue

                height_nn = dem_gk5_interpolate(br_gk5, bh_gk5)
                slope = slope_gk5_interpolate(br_gk5, bh_gk5)

                #print("bcol:", bcol, "crow/col:", (crow, ccol), "soil_id:", soil_id, "height_nn:", height_nn, "slope:", slope, "seed_harvest_cs:", seed_harvest_cs)

                #with open("dump-" + str(c) + ".json", "w") as jdf:
                #    json.dump({"id": (str(resolution) \
                #        + "|" + str(vrow) + "|" + str(vcol) \
                #        + "|" + str(crow) + "|" + str(ccol) \
                #        + "|" + str(soil_id) \
                #        + "|" + str(uj_id)), "sowing": worksteps[0], "harvest": worksteps[1]}, jdf, indent=2)
                #    c += 1

                # create the cutting worksteps according to the available grid data
                env_template["cropRotation"] = copy.deepcopy(crop_rotation_template)
                for year in range(cutting_year_start, cutting_year_end+1):
                    if year in mow_grids["no_of_cuts"]:
                        no_of_cuts = mow_grids["no_of_cuts"][year][(brow, bcol)]
                        if no_of_cuts > 0:
                            cm = copy.deepcopy(cutting_cm_template)
                            for cut_no in range(0, no_of_cuts + 1):
                                if cut_no < len(mow_grids["doys"][year]):
                                    rel_year = "000{cut_no}".format(cut_no=cut_no)
                                    cm["worksteps"][0]["date"] = cm["worksteps"][0]["date"].replace("0000", rel_year)
                                    cm["worksteps"][-1]["date"] = cm["worksteps"][0]["date"].replace("0000", rel_year)
                                    doy = mow_grids["doys"][year][cut_no][(brow, bcol)]
                                    d = date(year, 1, 1) + timedelta(days=int(doy)-1)
                                    cutting_ws = copy.deepcopy(cutting_ws_template)
                                    cutting_ws["date"] = "{rel_year}-{month:02d}-{day:02d}".format(rel_year=rel_year, month=d.month, day=d.day)
                                    cm["worksteps"].insert(-1, cutting_ws)
                            env_template["cropRotation"].append(cm)

                env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup["LeafExtensionModifier"]

                # set soil-profile
                soil_profile = monica_io3.find_and_replace_references(sp_json, sp_json)["result"]
                    
                #print("soil:", soil_profile)

                env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                # setting groundwater level
                if setup["groundwater-level"]:
                    groundwaterlevel = 20
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_in_groundwater", False):
                            groundwaterlevel = layer_depth
                            #print("setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m")
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [max(0, groundwaterlevel - 0.2) , "m"]
                    env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [groundwaterlevel + 0.2, "m"]
                    
                # setting impenetrable layer
                if setup["impenetrable-layer"]:
                    impenetrable_layer_depth = Mrunlib.get_value(env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_impenetrable", False):
                            impenetrable_layer_depth = layer_depth
                            #print("setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m")
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth, "m"]
                    env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

                if setup["elevation"]:
                    env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

                if setup["slope"]:
                    env_template["params"]["siteParameters"]["slope"] = slope / 100.0

                if setup["latitude"]:
                    clat, _ = cdict[(crow, ccol)]
                    env_template["params"]["siteParameters"]["Latitude"] = clat

                if setup["CO2"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(setup["CO2"])

                if setup["O3"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(setup["O3"])

                env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup["fertilization"]
                env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

                env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
                env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup["WaterDeficitResponseOn"]
                env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup["EmergenceMoistureControlOn"]
                env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup["EmergenceFloodingControlOn"]

                env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
                
                subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(climate_data=climate_data, 
                                                                  climate_model_folder=(climate_model + "/" if climate_model else ""),
                                                                  climate_scenario_folder=(climate_scenario + "/" if climate_scenario else ""),
                                                                  climate_region=climate_region,
                                                                  crow=str(crow),
                                                                  ccol=str(ccol))
                # subpath_to_csv = climate_data + "/csvs/" \
                # + (climate_model + "/" if climate_model else "") \
                # + (climate_scenario + "/" if climate_scenario else "") \
                # + climate_region + "/row-" + str(crow) + "/col-" + str(ccol) + ".csv"
                env_template["pathToClimateCSV"] = paths["monica-path-to-climate-dir"] + subpath_to_csv
                #print(env_template["pathToClimateCSV"])
                if DEBUG_WRITE_CLIMATE :
                    listOfClimateFiles.add(subpath_to_csv)

                env_template["customId"] = {
                    "setup_id": setup_id,
                    "brow": brow, "bcol": bcol,
                    "crow": int(crow), "ccol": int(ccol),
                    "soil_id": soil_id
                }

                if not DEBUG_DONOT_SEND :
                    socket.send_json(env_template)
                    print("sent env ", sent_env_count, " customId: ", env_template["customId"], flush=True)

                sent_env_count += 1

                # write debug output, as json file
                if DEBUG_WRITE:
                    if not os.path.exists(DEBUG_WRITE_FOLDER):
                        os.makedirs(DEBUG_WRITE_FOLDER)
                    if sent_env_count < DEBUG_ROWS  :

                        path_to_debug_file = DEBUG_WRITE_FOLDER + "/row_" + str(sent_env_count-1) + "_" + str(setup_id) + ".json" 

                        if not os.path.isfile(path_to_debug_file):
                            with open(path_to_debug_file, "w") as _ :
                                _.write(json.dumps(env_template))
                        else:
                            print("WARNING: Row ", (sent_env_count-1), " already exists", flush=True)
            #print("unknown_soil_ids:", unknown_soil_ids)

            #print("crows/cols:", crows_cols)
        stop_setup_time = time.clock()
        print("Setup ", (sent_env_count-1), " envs took ", (stop_setup_time - start_setup_time), " seconds", flush=True)

    stop_time = time.clock()

    # write summary of used json files
    if DEBUG_WRITE_CLIMATE:
        if not os.path.exists(DEBUG_WRITE_FOLDER):
            os.makedirs(DEBUG_WRITE_FOLDER)

        path_to_climate_summary = DEBUG_WRITE_FOLDER + "/climate_file_list" + ".csv"
        with open(path_to_climate_summary, "w") as _:
            _.write('\n'.join(listOfClimateFiles))

    try:
        print("sending ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds", flush=True)
        #print("ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
        print("exiting run_producer()", flush=True)
    except Exception:
        raise

if __name__ == "__main__":
    run_producer()