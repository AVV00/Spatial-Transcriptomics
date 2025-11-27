#!/usr/bin/env python
# coding=utf-8

#### import Moduules
import sys, os
import numpy as np
import csv
import gzip
import itertools
import tkinter
import zipfile
import re
from tkinter import messagebox

import matplotlib.pyplot as plt

import xopen
from isal import igzip, isal_zlib

import multiprocessing
from multiprocessing import Pool, cpu_count


###### Usage
usage = '''
     version %s
     Usage: %s <fqFile><outputPath><imageHeight><imageWidth> >STDOUT
''' % (PROG_VERSION,  os.path.basename(sys.argv[0]))

######  global variable  #####
p = re.compile(r'R\d{3}:C\d{3}')

### 解决tkinter 两个弹框的问题
# root = tkinter.Tk()
# root.withdraw()

def transCorrd(pos_x, pos_y, FovR, FovC, imageHeight, imageWidth):
    '''
    transfer the coord
    : pos_x - The original x coordinates
    : pos_y - The original y coordinates
    : FOVR - The row number of fov
    : FOVC - The col number of fov
    '''
    new_x = FovC * imageHeight - pos_y
    new_y = pos_x + (FovR - 1) * imageWidth
    return (round(new_x,2), round(new_y,2))
    #return (new_x, new_y)

def coordTransfer(fqFile, output, imageHeight, imageWidth, ratioHeight, ratioWidth, ratioHeightStart, ratioWidthStart, show = False):
    '''
    Convert Fov coordinates to chip coordinates
    '''
    outPutPath = os.path.join(output, "newCoord_" + os.path.splitext(os.path.basename(fqFile))[0] + ".gz" ) ## output Fastq File
    outFastq =  xopen.xopen(outPutPath, "wb", threads=0, compresslevel=3)

    if ratioHeightStart > ratioHeight:
        print("ratioHeightStart : %f is bigger than  ratioHeight : %f" % (ratioHeightStart, ratioHeight))
        return 0

    if ratioWidthStart > ratioWidth:
        print("ratioWidthStart : %f is bigger than  ratioWidth : %f" % (ratioWidthStart, ratioWidth))
        return 0

    ## Get central area information
    if ratioHeightStart < 0:
        ratioHeightStart = ratioHeight / 2

    if ratioWidthStart < 0:
        ratioWidthStart = ratioWidth / 2
    
    x_range = [round(imageWidth * ratioWidthStart), round(imageWidth * (1 - ratioWidth + ratioWidthStart))]
    y_range = [round(imageHeight * ratioHeightStart), round(imageHeight * (1 - ratioHeight + ratioHeightStart))]

    imageWidth_new = imageWidth * (1 - ratioWidth)  
    imageHeight_new = imageHeight * (1 - ratioHeight)
    
    if show:
        corrdRecord = []

    pre_row = 0
    with xopen.xopen(fqFile, threads=0) as pf:
        for line in pf:
            ### read infor for xopen.xopen
            title = line[:-1]
            read_seq = pf.readline()[:-1]
            Links = pf.readline()
            Q_value = pf.readline()[:-1]

            pres = p.findall(title) ## get the fov num of this read

            min_x = x_range[0]
            min_y = y_range[0]
            if len(pres) != 0:
                fov = pres[0]

                splitSet = title.split(" ")[-2].split(":")
                pos_y = float(splitSet[-1])
                pos_x = float(splitSet[-2])

                '''if int(fov[1:4]) == 1:
                    min_x = 0
                else:
                    min_x = x_range[0]

                if int(fov[5:8]) == 1:
                    min_y = 0
                else:
                    min_y = y_range[0]'''

                if (min_x < pos_x <= x_range[1]) & (min_y < pos_y <= y_range[1]): ##  Judge whether the coordinate is in the center of the image
                    
                    pos_x -= min_x
                    pos_y -= min_y
                    if int(fov[1:4]) != pre_row:
                        pre_row = int(fov[1:4])
                    new_pos = transCorrd(pos_x, pos_y, int(fov[1:4]), int(fov[6:9]), imageHeight_new, imageWidth_new)
                    if show:
                        corrdRecord.append(np.array(new_pos))

                    splitSet[-1] = str(new_pos[1])
                    splitSet[-2] = str(new_pos[0])
                    new_title = "_".join(splitSet)

                    outFastq.write((new_title+ '\n' + read_seq + '\n' + Links + Q_value + '\n').encode(encoding="utf-8"))
    
    outFastq.close()

    if show:
        corrdRecord = np.array(corrdRecord)
        corrdRecord = np.round(corrdRecord / 100).astype(int)
        corrdRecord = corrdRecord.T

        max_x = np.max(corrdRecord[0]) + 1
        max_y = np.max(corrdRecord[1]) + 1
        print(max_x, max_y)

        corrdRecord = corrdRecord.T

        heatmap = np.zeros((max_y, max_x))

        for readNum in range(len(corrdRecord)):
            coordTemp = corrdRecord[readNum]
            heatmap[coordTemp[1]][coordTemp[0]] += 1
        plt.matshow(heatmap)
        plt.savefig(os.path.join(output, 'heatMap.tif'))
        #cv2.imwrite("bins.tif", np.array(heatmap).astype(np.uint16))

def main():
    ######################### Phrase parameters #########################
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage)
    ArgParser.add_argument("--version", action="version", version=PROG_VERSION)
    ArgParser.add_argument("-r1", "--ratioWidth", action="store", dest="ratioWidth", type=float, default=0.0, metavar="FLOAT", help="Overlap ratio of width. [%(default)s]")
    ArgParser.add_argument("-r1s", "--ratioWidthStart", action="store", dest="ratioWidthStart", type=float, default=-1.0, metavar="FLOAT", help="start Overlap ratio of width. [%(default)s]")
    ArgParser.add_argument("-r2", "--ratioHeight", action="store", dest="ratioHeight", type=float, default=0.0, metavar="FLOAT", help="Overlap ratio of heigh [%(default)s]")
    ArgParser.add_argument("-r2s", "--ratioHeightStart", action="store", dest="ratioHeightStart", type=float, default=-1.0, metavar="FLOAT", help="start Overlap ratio of heigh [%(default)s]")
    ArgParser.add_argument("-s", "--showHeatmap", action="store_true", dest="showHeatmap", default=True, help="Display heatmap. [%(default)s]")

    (para, args) = ArgParser.parse_known_args()

    if len(args) != 4:
        ArgParser.print_help()
        print(sys.stderr, "\nERROR: The parameters number is not correct!")
        sys.exit(1)
    else:
        (fqFile, output, imageHight, imageWidth) = args

    
    coordTransfer(fqFile, output, int(imageHight), int(imageWidth), float(para.ratioHeight), float(para.ratioWidth), float(para.ratioHeightStart), float(para.ratioWidthStart), bool(para.showHeatmap))

if __name__ == "__main__":
    main()
    

