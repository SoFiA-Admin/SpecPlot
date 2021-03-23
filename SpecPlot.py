#!/usr/bin/env python
### ____________________________________________________________________ ###
###                                                                      ###
### SpecPlot.py - SpecPlot 0.1 - Plotting tool for SoFiA output spectra  ###
### Copyright (C) 2021 Tobias Westmeier                                  ###
### ____________________________________________________________________ ###
###                                                                      ###
### Address:  Tobias Westmeier                                           ###
###           ICRAR M468                                                 ###
###           The University of Western Australia                        ###
###           35 Stirling Highway                                        ###
###           Crawley WA 6009                                            ###
###           Australia                                                  ###
###                                                                      ###
### E-mail:   tobias.westmeier [at] uwa.edu.au                           ###
### ____________________________________________________________________ ###
###                                                                      ###
### This program is free software: you can redistribute it and/or modify ###
### it under the terms of the GNU General Public License as published by ###
### the Free Software Foundation, either version 3 of the License, or    ###
### (at your option) any later version.                                  ###
###                                                                      ###
### This program is distributed in the hope that it will be useful,      ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of       ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ###
### GNU General Public License for more details.                         ###
###                                                                      ###
### You should have received a copy of the GNU General Public License    ###
### along with this program. If not, see http://www.gnu.org/licenses/.   ###
### ____________________________________________________________________ ###
###                                                                      ###


# ---------------------------------------------
# Import Python modules
# ---------------------------------------------
import sys
import os
import math
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


# ---------------------------------------------
# Global settings
# ---------------------------------------------
# These can be adjusted by the user to modify
# the appearance of the plots. Global settings
# for matplotlib are described here:
#
# https://matplotlib.org/users/customizing.html
# ---------------------------------------------
plt.rcParams["figure.figsize"]         =  (9.0,6.0);
mpl.rcParams["axes.linewidth"]         =  0.75;
mpl.rcParams["axes.labelweight"]       =  "light";
mpl.rcParams["xtick.top"]              =  True;
mpl.rcParams["xtick.bottom"]           =  True;
mpl.rcParams["ytick.left"]             =  True;
mpl.rcParams["ytick.right"]            =  True;
mpl.rcParams["xtick.direction"]        =  "in";
mpl.rcParams["ytick.direction"]        =  "in";
mpl.rcParams["xtick.major.width"]      =  0.5;
mpl.rcParams["ytick.major.width"]      =  0.5;
mpl.rcParams["grid.linewidth"]         =  0.25;
mpl.rcParams["grid.linestyle"]         =  "solid";
mpl.rcParams["grid.color"]             =  "#E0E0E0";
mpl.rcParams["grid.alpha"]             =  1.0;
mpl.rcParams["lines.linewidth"]        =  0.75;
mpl.rcParams["lines.solid_joinstyle"]  =  "miter";
mpl.rcParams["lines.solid_capstyle"]   =  "projecting";
mpl.rcParams["patch.linewidth"]        =  0.75;
mpl.rcParams["font.family"]            =  "sans-serif";
mpl.rcParams["font.size"]              =  11;
mpl.rcParams["font.weight"]            =  "light";
color_spectrum                         =  "#2040A0";
color_uncertainty                      =  "#E0E0E0";


# ---------------------------------------------
# A few natural constants
# ---------------------------------------------
speed_of_light = 299792.458;  #km/s
restfreq = 1420.4057517667;   # MHz


# ---------------------------------------------
# Check command line arguments
# ---------------------------------------------
if len(sys.argv) < 3 or len(sys.argv) > 10:
	sys.stderr.write("\nUsage: SpecPlot.py <input> <output> [<velo>] [<rms> <beam>]\n       [x_min x_max y_min y_max]\n\n");
	sys.stderr.write("       <input>   Name of the input spectrum created by SoFiA.\n");
	sys.stderr.write("       <output>  Name of the output plot to be written to disc. Must\n");
	sys.stderr.write("                 end in '.eps' or '.pdf'.\n");
	sys.stderr.write("       <velo>    Frequency-to-velocity conversion method. Can be 'cz'\n");
	sys.stderr.write("                 or 'rest' for recession velocity or source rest-frame\n");
	sys.stderr.write("                 velocity, respectively.\n");
	sys.stderr.write("       <rms>     Local rms noise level in the data cube from which the\n");
	sys.stderr.write("                 spectrum was extracted (in the same units as the flux\n");
	sys.stderr.write("                 values in the original spectrum).\n");
	sys.stderr.write("       <beam>    Solid angle of the telescope beam in units of pixels.\n");
	sys.stderr.write("       <x_min>   Plotting range lower limit in x.\n");
	sys.stderr.write("       <x_max>   Plotting range upper limit in x.\n");
	sys.stderr.write("       <y_min>   Plotting range lower limit in y.\n");
	sys.stderr.write("       <y_max>   Plotting range upper limit in y.\n\n");
	
	sys.stderr.write("       <input> and <output> specify the names of the input spectrum and\n");
	sys.stderr.write("       the output plot, respectively. Both parameters are mandatory.\n");
	sys.stderr.write("       The output file type will be determined from the file ending and\n");
	sys.stderr.write("       can either be PDF (file ending '.pdf') or EPS (file ending '.eps').\n\n");
	
	sys.stderr.write("       If <velo> is specified, then the frequency axis of the spectrum\n");
	sys.stderr.write("       will be automatically converted to velocity. Possible values\n");
	sys.stderr.write("       for <velo> are 'cz' or 'rest' for recession velocity or source\n");
	sys.stderr.write("       rest-frame velocity, respectively. Alternatively, if set to\n");
	sys.stderr.write("       'freq', then no conversion will be done. If 'rest' is selected,\n");
	sys.stderr.write("       then the centroid of the spectrum will define v = 0. Note that\n");
	sys.stderr.write("       velocity conversion will only work correctly for the 21-cm line\n");
	sys.stderr.write("       of neutral atomic hydrogen, as the restfrequency is hardcoded!\n\n");
	
	sys.stderr.write("       <rms> and <beam> must be specified for any beam correction to\n");
	sys.stderr.write("       take place and for statistical uncertainties to be plotted as\n");
	sys.stderr.write("       well. For a Gaussian beam, the solid angle can be calculated\n");
	sys.stderr.write("       as (pi / 4) * a * b / ln(2) where a and b are the major and\n");
	sys.stderr.write("       minor axis in units of pixels.\n\n");
	
	sys.stderr.write("       The <x_min>, <x_max>, <y_min> and <y_max> parameters can be\n");
	sys.stderr.write("       used to specify the plotting range in the actual units in which\n");
	sys.stderr.write("       the data will appear in the plot. It might be useful to first\n");
	sys.stderr.write("       leave the range unspecified to check the correct output range\n");
	sys.stderr.write("       and units before setting the range limits.\n\n");
	sys.exit(0);

velotype = "freq";
rms      = 0.0;
beam     = 1.0;
x_min    = 0.0;
x_max    = 0.0;
y_min    = 0.0;
y_max    = 0.0;

input_file  = str(sys.argv[1]);
output_file = str(sys.argv[2]);

if len(sys.argv) == 4:
	velotype = str(sys.argv[3]);
elif len(sys.argv) == 5:
	rms      = float(sys.argv[3]);
	beam     = float(sys.argv[4]);
elif len(sys.argv) == 6:
	velotype = str(sys.argv[3]);
	rms      = float(sys.argv[4]);
	beam     = float(sys.argv[5]);
elif len(sys.argv) == 7:
	x_min    = float(sys.argv[3]);
	x_max    = float(sys.argv[4]);
	y_min    = float(sys.argv[5]);
	y_max    = float(sys.argv[6]);
elif len(sys.argv) == 8:
	velotype = str(sys.argv[3]);
	x_min    = float(sys.argv[4]);
	x_max    = float(sys.argv[5]);
	y_min    = float(sys.argv[6]);
	y_max    = float(sys.argv[7]);
elif len(sys.argv) == 9:
	rms      = float(sys.argv[3]);
	beam     = float(sys.argv[4]);
	x_min    = float(sys.argv[5]);
	x_max    = float(sys.argv[6]);
	y_min    = float(sys.argv[7]);
	y_max    = float(sys.argv[8]);
elif len(sys.argv) == 10:
	velotype = str(sys.argv[3]);
	rms      = float(sys.argv[4]);
	beam     = float(sys.argv[5]);
	x_min    = float(sys.argv[6]);
	x_max    = float(sys.argv[7]);
	y_min    = float(sys.argv[8]);
	y_max    = float(sys.argv[9]);

if velotype != "freq" and velotype != "cz" and velotype != "rest":
	sys.stderr.write("ERROR: Unknown velocity type; must be 'freq', 'cz' or 'rest'.\n");
	sys.exit(1);


# ---------------------------------------------
# Establish output format
# ---------------------------------------------
output_format = (output_file[-3:]).lower();
if output_format != "eps" and output_format != "pdf":
	sys.stderr.write("ERROR: Unknown output format; must be 'eps' or 'pdf'.\n");
	sys.exit(1);


# ---------------------------------------------
# Try to read input file
# ---------------------------------------------
column_names = [];
column_units = [];
try:
	# Read data table
	with open(input_file, "r") as fp:
		spectrum = np.loadtxt(fp, dtype="float", comments="#", unpack=True);
	
	# Read header info
	with open(input_file, "r") as fp:
		lines = fp.readlines();
		for line in lines:
			if line[0:10] == "#  Channel":
				column_names.append(line[1:10].strip());
				column_names.append(line[10:28].strip());
				column_names.append(line[28:46].strip());
				column_names.append(line[46:56].strip());
			elif line[0:10] == "#        -":
				column_units.append(line[1:10].strip());
				column_units.append(line[10:28].strip());
				column_units.append(line[28:46].strip());
				column_units.append(line[46:56].strip());
				break;
except:
	sys.stderr.write("ERROR: Failed to read spectrum:\n  {:s}\n".format(input_file));
	sys.exit(1)


# ---------------------------------------------
# Extract data columns
# ---------------------------------------------
xaxis = spectrum[1];
data  = np.nan_to_num(spectrum[2]);
error = rms * np.sqrt(spectrum[3] / beam);


# ---------------------------------------------
# Do some unit conversions for convenience
# ---------------------------------------------
if column_units[1] == "m/s":
	column_units[1] = "km/s";
	xaxis /= 1e+3;
elif column_units[1] == "Hz":
	column_units[1] = "MHz";
	xaxis /= 1e+6;

if column_units[2] == "Jy/beam" and rms > 0.0 and beam > 0.0:
	column_units[2] = "mJy";
	data *= 1e+3 / beam;
	error *= 1e+3;
elif column_units[2] == "Jy":
	column_units[2] = "mJy";
	data *= 1e+3;
	error *= 1e+3;


# ---------------------------------------------
# Convert frequency to recession velocity or
# source rest frame velocity
# ---------------------------------------------
if column_units[1] == "MHz":
	if velotype == "rest":
		sys.stderr.write("WARNING: Converting frequency to source rest-frame velocity.\n");
		centroid = 0.0;
		counter = 0.0;
		for i in range(len(xaxis)):
			centroid += xaxis[i] * data[i];
			counter += data[i];
		xaxis = speed_of_light * (1.0 - xaxis * counter / centroid);
		column_names[1] = "Rest-frame velocity"
		column_units[1] = "km/s"
	elif velotype == "cz":
		sys.stderr.write("WARNING: Converting frequency to recession velocity (cz).\n");
		xaxis = speed_of_light * ((restfreq / xaxis) - 1.0);
		column_names[1] = "Recession velocity"
		column_units[1] = "km/s"


# ---------------------------------------------
# Set scale
# ---------------------------------------------
upper = data + error;
lower = data - error;
if x_min >= x_max or y_min >= y_max:
	x_min = np.nanmin(xaxis);
	x_max = np.nanmax(xaxis);
	y_min = np.nanmin(lower);
	y_max = np.nanmax(upper);
	y_min -= 0.1 * (y_max - y_min);
	y_max += 0.1 * (y_max - y_min);


# ---------------------------------------------
# Ensure that baseline extends all the way to
# edge if possible
# ---------------------------------------------
if data[0] == 0.0 and data[-1] == 0.0:
	if xaxis[0] < xaxis[-1]:
		xaxis[0]  = min(xaxis[0],  x_min);
		xaxis[-1] = max(xaxis[-1], x_max);
	else:
		xaxis[0]  = max(xaxis[0],  x_max);
		xaxis[-1] = min(xaxis[-1], x_min);


# ---------------------------------------------
# Create plot
# ---------------------------------------------
if rms > 0.0: plt.fill_between(xaxis, lower, upper, color=color_uncertainty, step="mid", joinstyle="miter");
plt.step(xaxis, data, where="mid", color=color_spectrum);
plt.xlabel("{:s} ({:s})".format(column_names[1], column_units[1]));
plt.ylabel("{:s} ({:s})".format(column_names[2], column_units[2]));
plt.grid(True);
plt.xlim([x_min, x_max]);
plt.ylim([y_min, y_max]);
plt.savefig(output_file, orientation="portrait", format=output_format, papertype="a4", bbox_inches="tight");
