import pandas as pd
import numpy as np
import os
import cv2
import datetime

from bs4 import BeautifulSoup

# inner: is the metadata inside the acq_folder or outside?
def collect_metadata(xml_path, inner=False):    

    # Read xml and prepare for the soup
    with open(xml_path) as fp:
        soup = BeautifulSoup(fp, "xml")

    # Check for merged?
    # Collect info about the tilescan
    tilescan_info = soup.find_all("Tile")
    # Variables to keep the information
    xix_lst = np.zeros(len(tilescan_info), dtype=np.int32)
    yix_lst = np.zeros_like(xix_lst)
    xpos_lst = np.zeros_like(xix_lst, dtype=np.double)
    ypos_lst = np.zeros_like(yix_lst, dtype=np.double)
    # Run through each tile and save the position
    for tile_idx in range(len(tilescan_info)):
        tile = tilescan_info[tile_idx]
        xix_lst[tile_idx] = tile.get("FieldX")
        yix_lst[tile_idx] = tile.get("FieldY")
        xpos_lst[tile_idx] = tile.get("PosX")
        ypos_lst[tile_idx] = tile.get("PosY")

    xix_unique_ar = np.unique(xix_lst)
    yix_unique_ar = np.unique(yix_lst)
    tiles = {"xix_lst": xix_lst,
            "yix_lst": yix_lst,
            "xix_unique_ar": xix_unique_ar,
            "yix_unique_ar": yix_unique_ar,
            "xpos_lst": xpos_lst*1e6,
            "ypos_lst": ypos_lst*1e6,
            "tile_xcnt": len(xix_unique_ar),
            "tile_ycnt": len(yix_unique_ar)}

    # Collect dimenions: x,y,z,t,stage (stage gives info about tilescan)
    dimension_desc = soup.find_all("DimensionDescription")
    dimensions = {}
    for desc in dimension_desc:
        dimid = desc.get("DimID")
        unit = desc.get("Unit").replace("\xb5", "u")[-2:] # Replace greek letter of \mu with letter u
        length = desc.get("Length")
        counts = np.int32(desc.get("NumberOfElements"))
        voxel = desc.get("Voxel")
        # convert all length values to um
        if len(unit) > 0:   # implies a numerical value
            length = np.double(length)
            voxel = np.double(voxel)
            if unit == "mm":
                length *= 1000
                voxel *= 1000
                unit = "um"
        dimensions[dimid] = {"Length": length,
                             "NumberOfElements": counts,
                             "Unit": unit,
                             "Voxel": voxel
                            }
    
    # Datetime object from the StartTime
    datetime_str = soup.find("StartTime").text
    datetime_obj = datetime.datetime.strptime(datetime_str, "%m/%d/%Y %I:%M:%S %p.%f")
    return {"dimensions": dimensions,
            "tiles": tiles,
            "xml_path": xml_path,
            "start_time": datetime_obj}
