import pandas as pd
import numpy as np
import os
import cv2
import datetime

from bs4 import BeautifulSoup

# exp_name = "EQ59_Gly_02022021"
# acq_name = "TileScan_Tp1-5"

# exp_folder = f"D:/Tolga/Colony Images/{exp_name}/"
# base_folder = exp_folder + f"{exp_name}_{acq_name}/"
# metadata_folder = base_folder + "MetaData/"

# xml_path = metadata_folder + f"{exp_name}_{acq_name}_Properties.xml"

# inner: is the metadata inside the acq_folder or outside?
def collect_metadata(exp_folder, exp_name, acq_name, inner=False):
    tilescan_folder = os.path.join(exp_folder, exp_name)
    if inner:
        acq_folder = os.path.join(tilescan_folder, f"{acq_name}")
        metadata_folder = os.path.join(tilescan_folder, f"{acq_name}", "MetaData")
        xml_path = os.path.join(metadata_folder, f"{acq_name}_Properties.xml")
    else:
        acq_folder = os.path.join(tilescan_folder, f"{exp_name}_{acq_name}")
        metadata_folder = os.path.join(tilescan_folder, "MetaData")
        xml_path = os.path.join(metadata_folder, f"{acq_name}_Properties.xml")
    

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
            "exp_folder": exp_folder,
            "acq_folder": acq_folder,
            "exp_name": exp_name,
            "acq_name": acq_name,
            "xml_path": xml_path,
            "start_time": datetime_obj}
            


# Block is defined by the x and y values of each edge
class Block:
    # __init__(self, left_top, right_bottom, width, height, xvoxel, yvoxel):
    def __init__(self, left_top, right_bottom, width, height, xvoxel, yvoxel):
        """
            @param left_top: The location of the left_top point
            @param right_bottom: The location of the right bottom
            @param width: Number of pixels in x
            @param height: Number of pixels in y
            @param xvoxel: width of each pixel in physical length
            @param yvoxel: height of each pixel in physical length
            
            left_top and right_bottom are calculated from the xml.properties file inside the Metadata folder
            The edges self.left, self.right, self.top and self.bottom defines the corresponding edges
            
            self.xidx and self.yidx are the indices that would correspond to the image
            
            Purpose of this class:
            The tilescan images have overlapping grids, so to merge all the grid-based tilescan images, 
            we need to redefine the block in the positions of each grid block
            
            How to use:
            First, we initialize each block with the overlapping coordinates. Then, consequentially, fix the overlapping
            coordinates by using the function overlap_fix(block). The parameter block here is an object of Block, which is near.
            The overlap_fix first checks where the overlap occurs with the block, then redefines the overlapping edge
            
            Note: The overlap_fix function consequentially assume overlapping blocks have the same non-overlapping position.
                Make sure to check this with the plot below.
        """
        # Edges
        self.left = left_top[0]
        self.right = right_bottom[0]
        self.top = left_top[1]
        self.bottom = right_bottom[1]
        
        # x,y pixel count
        self.width = width
        self.height = height
        # x and y pixel size in physical length
        self.xvoxel = xvoxel
        self.yvoxel = yvoxel
        # x,y pixel indices
        self.xidx = np.arange(self.width)
        self.yidx = np.arange(self.height)
        
    # def is_inside_h(self, y):
    #     # Checks whether the y-coordinate is inside the horizontal range defined by the block
    #     if self.bottom < y and self.top > y:
    #         return True
    #     else:
    #         return False
    # def is_inside_v(self, x):
    #     # Checks whether the x-coordinate is inside the vertical range defined by the block
    #     if self.left < x and self.right > x:
    #         return True
    #     else:
    #         return False
    # def is_inside(self, x, y):
    #     """
    #         There must be a check for the microscope parameters flipx and swap-xy
    #         The function returns true if the given location is inside the block
            
    #     """
    #     print("Comparing x- (%g,%g) -- %g" % (left,right, x))
    #     print("Comparing y- (%g,%g) -- %g" % (bottom, top, y))
    #     if self.left < x and self.right > x and self.bottom < y and self.top > y: # inside the block
    #         return True
    #     else:
    #         return False
    def remove_left(self, right, xsz): # Remove the overlapping pixels up to right
        # The length of the overlap
        length = right - self.left
        # print("x - Length = ", length)
        # self.left = right
        # Overlapping is removed
        self.xidx = np.arange(int(length/self.xvoxel), xsz)
        # print("First x-index = ", self.xidx[0])
    def remove_bottom(self, top, ysz): # Remove the overlapping pixels up to top
        # The length of the overlap
        length = top - self.bottom
        # print("y - Length = ", length)
        # selt.
        # Overlapping is removed
        self.yidx = np.arange(int(length/self.yvoxel), ysz)
        # print("First y-index = ", self.yidx[0])
        
    # def overlap_fix(self, block):
    #     # Overlap on the right --> Right of the self block is the left of the param block
    #     if self.is_inside_v(block.left):
    #         block.remove_left(self.right)
    #     if self.is_inside_h(block.bottom):
    #         block.remove_bottom(self.top)
        

def merge(metadata, channel, inner=False): # For now, assumes all the acquire is xyzt

    ###############################
    #### Construct the filepaths for each tilescan tif image
    ###############################
    # Construct the filepaths to collect the tilescan images based on the provided metadata
    znum = metadata["dimensions"]["Z"]["NumberOfElements"] # Number of z planes
    # Find the number of digits
    znum_digit = len(str(znum-1))
    # zstr for file path
    zstr_holder = f"z%0{znum_digit}d"

    snum = metadata["dimensions"]["Stage"]["NumberOfElements"] # Number of tilescan images
    # Find the number of digits
    snum_digit = len(str(snum-1))
    # sstr for file path
    sstr_holder = f"s%0{snum_digit}d"

    # Check if tstr exists
    tstr_exists = True
    try:
        tnum = metadata["dimensions"]["T"]["NumberOfElements"] # Number of tilescan images
        # Find the number of digits
    except KeyError:
        # tnum is not set
        tnum = 1
        tstr_exists = False

    tnum_digit = len(str(tnum-1))
    # sstr for file path
    tstr_holder = f"t%0{tnum_digit}d"        

    # folder and filenames
    exp_folder = metadata["exp_folder"]
    acq_folder = metadata["acq_folder"]
    exp_name = metadata["exp_name"]
    acq_name = metadata["acq_name"]

    # fname_list = []
    # for tix in range(tnum):
    #     for six in range(snum):
    #         for zix in range(znum):
    #             tstr = tstr_holder % (tix)
    #             sstr = sstr_holder % (six)
    #             zstr = zstr_holder % (zix)
    #             # Combine all params
    #             fname = f"{base_folder}TileScan/{exp_name}_{acq_name}_{tstr}_{sstr}_{zstr}.tif"
    #             fname_list.append(fname)
    ###########################################
    ###### Calculating the positions of each tile
    ###### Constructing the blocks
    ###########################################
    xix_unique_ar = metadata["tiles"]["xix_unique_ar"]
    yix_unique_ar = metadata["tiles"]["yix_unique_ar"]
    xpos_lst = metadata["tiles"]["xpos_lst"]
    ypos_lst = metadata["tiles"]["ypos_lst"]
    xix_lst = metadata["tiles"]["xix_lst"]
    yix_lst = metadata["tiles"]["yix_lst"]
    xvoxel = metadata["dimensions"]["X"]["Voxel"]
    yvoxel = metadata["dimensions"]["Y"]["Voxel"]
    xsz = metadata["dimensions"]["X"]["NumberOfElements"]
    ysz = metadata["dimensions"]["Y"]["NumberOfElements"]
    # Define all the blocks and plot the edges
    blocks_2d = np.empty((len(xix_unique_ar), len(yix_unique_ar)), dtype=Block)
    for six in range(snum):
        (x,y) = (xpos_lst[six], ypos_lst[six])   # convert to um
        # Define the square
        left = x-xvoxel*xsz/2
        right = x+xvoxel*xsz/2
        top = y+yvoxel*ysz/2
        bottom = y-yvoxel*ysz/2
        
        xix = xix_lst[six]
        yix = yix_lst[six]
        
        blocks_2d[xix, yix] = Block((left, top), (right, bottom), xsz, ysz, xvoxel, yvoxel)
    ##########################################
    ####### Remove overlap from neighbor blocks
    ##########################################
    for xix in xix_unique_ar[1:]:
        for yix in yix_unique_ar:
            # remove vertical overlap
            blocks_2d[xix, yix].remove_left(blocks_2d[xix-1,yix].right, xsz)
    for yix in yix_unique_ar[1:]:
        for xix in xix_unique_ar:
            # remove horizontal overlap
            blocks_2d[xix, yix].remove_bottom(blocks_2d[xix,yix-1].top, ysz)
    ################################
    ####### Merge using the overlap-removed blocks
    ################################
    merged_folder = os.path.join(exp_folder, "Merged/Channel_"+channel)
    # Let's test merging a single timepoint for each z sections
    # Run through all the tstr and zstr
    tstr = tstr_holder % (0)
    zstr = zstr_holder % (0)
    sstr = sstr_holder % (0)
    # fpath_test = os.path.join(acq_folder, f"{exp_name}_{acq_name}_{tstr}_{sstr}_{zstr}_ch00.tif")
    # img_test = cv2.imread(fpath_test, cv2.IMREAD_GRAYSCALE)
    # test_height, test_width = img_test.shape
    
    test_height = ysz
    test_width = xsz
    
    # For each timepoint
    for tix in range(tnum):
        # Only if tstr_holder exists
        tstr = tstr_holder % (tix)
        
        # For each z section
        for zix in range(znum):
            zstr = zstr_holder % (zix)
            # File path to merge into
            merged_acq_folder = os.path.join(merged_folder, f"{exp_name}_{acq_name}")
            if not os.path.isdir(merged_acq_folder):
                os.mkdir(merged_acq_folder)
            
            fpath_merged = os.path.join(merged_acq_folder, f"{exp_name}_{acq_name}_Merged_{tstr}_{zstr}.tif")

            # Create empty merged image
            # img_merged_bw = np.zeros((test_height*len(yix_unique_ar), test_width*len(xix_unique_ar)), dtype=img_test.dtype)
            img_merged_bw = np.zeros((test_height*len(yix_unique_ar), test_width*len(xix_unique_ar)), dtype=np.uint8)

            for six in range(snum):    # Run over the stage positions -> single merged tif file
                sstr = sstr_holder % (six)
                #Construct the filepath to merge
                #fpath = getFileName(fname_list[0], zstr_ix, zstr, tstr_ix, tstr, sstr_ix, sstr)

                if channel == "gfp":
                    if inner:
                        if tstr_exists:
                            fpath = os.path.join(acq_folder, f"{acq_name}_{tstr}_{sstr}_{zstr}_ch00.tif")
                        else:
                            fpath = os.path.join(acq_folder, f"{acq_name}_{sstr}_{zstr}_ch00.tif")
                    else:
                        if tstr_exists:
                            fpath = os.path.join(acq_folder, f"{exp_name}_{acq_name}_{tstr}_{sstr}_{zstr}_ch00.tif")
                        else:
                            fpath = os.path.join(acq_folder, f"{exp_name}_{acq_name}_{sstr}_{zstr}_ch00.tif")
                else:
                    if inner:
                        if tstr_exists:
                            fpath = os.path.join(acq_folder, f"{acq_name}_{tstr}_{sstr}_{zstr}_ch01.tif")
                        else:
                            fpath = os.path.join(acq_folder, f"{acq_name}_{sstr}_{zstr}_ch01.tif")
                    else:
                        if tstr_exists:
                            fpath = os.path.join(acq_folder, f"{exp_name}_{acq_name}_{tstr}_{sstr}_{zstr}_ch01.tif")
                        else:
                            fpath = os.path.join(acq_folder, f"{exp_name}_{acq_name}_{sstr}_{zstr}_ch01.tif")   

                if not os.path.isfile(fpath):
                    print(f"TileScan image file path not defined:\n\t {fpath}")
                
                # x,y indices
                xix = xix_lst[six]
                yix = yix_lst[six]
                block = blocks_2d[xix,yix]
                
                img = cv2.imread(fpath, cv2.IMREAD_GRAYSCALE)
                #height,width = img.shape
                width = len(block.xidx)
                height = len(block.yidx)
            
                # If the img is empty, initiate with all zeros
                if img is None:
                    img = np.zeros((test_width, test_height))
                
                # # Image pixel positions
                # xixar = np.arange(0,width) + width*xix
                # yixar = np.arange(0,height) + height*yix
                
                # SwapXY and FlipX
                flipped_img = np.transpose(np.flipud(img))
                #####################################################################
                #####################################################################
                #img_merged_bw[(yix*height):((yix+1)*height), (xix*width):((xix+1)*width)] = flipped_img[block.yidx[0]:(block.yidx[-1]+1), block.xidx[0]:(block.xidx[-1]+1)]

                if xix == 0 and yix == 0:
                    img_merged_bw[(yix*height):((yix+1)*height), (xix*width):((xix+1)*width)] = flipped_img[block.yidx[0]:(block.yidx[-1]+1), block.xidx[0]:(block.xidx[-1]+1)]
                elif xix == 0 and yix !=0:
                    img_merged_bw[ysz+((yix-1)*height):ysz+(yix*height), (xix*width):((xix+1)*width)] = flipped_img[block.yidx[0]:(block.yidx[-1]+1), block.xidx[0]:(block.xidx[-1]+1)]
                elif xix != 0 and yix ==0:
                    img_merged_bw[(yix*height):((yix+1)*height), xsz+((xix-1)*width):xsz+(xix*width)] = flipped_img[block.yidx[0]:(block.yidx[-1]+1), block.xidx[0]:(block.xidx[-1]+1)]
                elif xix != 0 and yix !=0:
                    img_merged_bw[ysz+((yix-1)*height):ysz+(yix*height), xsz+((xix-1)*width):xsz+(xix*width)] = flipped_img[block.yidx[0]:(block.yidx[-1]+1), block.xidx[0]:(block.xidx[-1]+1)]

            print(f"Image merged for zstr = {zstr} ({zix}/{znum}), tstr = {tstr} ({tix}/{tnum})",
                    end="\r", flush=True)

            # Write the merged image
            cv2.imwrite(fpath_merged, img_merged_bw)

def make_movie(metadata, channel, ROIx, ROIy, RR, dim_vid, ksz_um, marked=False, full_stage=False):
    # For a set number of t
    # dim_desc["tsz"] = 1
    # For a set number of z
    # dim_desc["zsz"] = 80
    xunit = metadata["dimensions"]["X"]["Unit"]
    zunit = metadata["dimensions"]["Z"]["Unit"]

    xvoxel = metadata["dimensions"]["X"]["Voxel"]
    zvoxel = metadata["dimensions"]["Z"]["Voxel"]
    
    # tsz = metadata["dimensions"]["T"]["NumberOfElements"]
    xsz = metadata["dimensions"]["X"]["NumberOfElements"]
    ysz = metadata["dimensions"]["Y"]["NumberOfElements"]
    zsz = metadata["dimensions"]["Z"]["NumberOfElements"]

    w = metadata["tiles"]["xix_unique_ar"].size*xsz
    h = metadata["tiles"]["yix_unique_ar"].size*ysz

    ksz_um = 10 # Kernel size for filters in um

    ksz_px = ksz_um/xvoxel
    # Make kernel size an odd integer
    ksz_px = int(ksz_px/2)*2 + 1
    print(f"width: {w} - height: {h}")
    print(f"kernel um: {ksz_um} - kernel px: {ksz_px}")

    # xix_unique_ar = metadata["tiles"]["xix_unique_ar"]
    # xsz = metadata["dimensions"]["X"]["NumberOfElements"]

    exp_name = metadata["exp_name"]
    acq_name = metadata["acq_name"]
    exp_folder = metadata["exp_folder"]

    if marked:
        marked_str = "Marked"
    else:
        marked_str = "Unmarked"

    if full_stage:
        full_stage_str = "Full"
    else:
        full_stage_str = "ROI"

    video_folder = os.path.join(exp_folder, f"Videos_{full_stage_str}_{marked_str}/Channel_"+channel)

    merged_folder = os.path.join(exp_folder, "Merged/Channel_"+channel)

    # # Find the number of digits
    # tnum_digit = len(str(tsz))
    # # tstr for file path
    # tstr_holder = f"t%0{tnum_digit}d"

    # Find the number of digits
    znum_digit = len(str(zsz-1))
    # sstr for file path
    zstr_holder = f"z%0{znum_digit}d"
    
    fpath_test = os.path.join(merged_folder, f"{exp_name}_{acq_name}", 
                f"{exp_name}_{acq_name}_Merged_t0_{zstr_holder % 0}.tif")
    img_test = cv2.imread(fpath_test, cv2.IMREAD_GRAYSCALE)
    if full_stage:
        test_height = img_test.shape[0]
    else:
        test_height = 2*RR

    scale = dim_vid[0]/test_height

    video_acq_folder = os.path.join(video_folder, f"{exp_name}_{acq_name}")
    if not os.path.isdir(video_acq_folder):
        os.mkdir(video_acq_folder)
    video_path = os.path.join(video_acq_folder, f"{exp_name}_{acq_name}.avi")
    out = cv2.VideoWriter(video_path, cv2.VideoWriter_fourcc(*'DIVX'), 4, dim_vid)
    # img_zsum = np.zeros(zsz)
    for zix in range(zsz):
        zstr = zstr_holder % (zix)
        merged_path = os.path.join(merged_folder, f"{exp_name}_{acq_name}", f"{exp_name}_{acq_name}_Merged_t0_{zstr}.tif")

        print(f"Loading file = {merged_path}", end="\r", flush=True)
        
        img = cv2.imread(merged_path)
        # Choose ROI image
        if full_stage:
            roi_img = img
        else:
            roi_img = img[ROIy[0]:ROIy[1], ROIx[0]:ROIx[1]]
        img_blur = cv2.GaussianBlur(roi_img, (ksz_px,ksz_px), 0)

        # Add contrast for better visualization.
        alpha = 1
        beta = 0
        img_contrast = np.uint8(np.clip(alpha*img_blur + beta, 0, 255))
        # Resize for video output
        img_resized = cv2.resize(img_contrast, dim_vid)

        if marked:
            # Scale bar
            scalebar_length = 500   # in um
            bar_start_coor = (50,440)
            bar_end_coor = (int(50 + scalebar_length/xvoxel*scale),450)
            bar_thickness = -1
            bar_color = (255,255,255)
            cv2.rectangle(img_resized, bar_start_coor, bar_end_coor, bar_color, bar_thickness)
            
            text_x = int(bar_start_coor[0])
            text_y = int((bar_start_coor[1]+bar_end_coor[1])/2 + 30)
            cv2.putText(img_resized,"%d %s" % (scalebar_length, xunit), (text_x,text_y), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,255,255), 2)
            
            zval_x = 260
            zval_y = 50
            z_str = "z=%.3f%s" % (zix*zvoxel, zunit)
            cv2.putText(img_resized, z_str, (zval_x, zval_y), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,255,255), 2)
            
            # tval_x = 0
            # tval_y = 20
            # t_str_show = "t=20h %d min" % (27)
            # cv2.putText(img_jpg, t_str_show, (tval_x, tval_y), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,255,255), 2)
            
        out.write(img_resized)
        # img_zsum[zix] = sum(sum(sum(img_blur/255.)))/3
        
    out.release()