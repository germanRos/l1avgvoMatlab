﻿l1avgvoMatlab
=============

Fast and Robust l1-averaging-based Pose Estimation for Driving Scenarios

 Copyright 2013. All rights reserved.
 Author: German Ros (gros@cvc.uab.es)
         Advanced Driver Assistance Systems (ADAS)
         Computer Vision Center (CVC)
	  Universitat Autònoma de Barcelona (UAB)

 This is free software; you can redistribute it and/or modify it under the
 terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or any later version.

 This software is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS     
 FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this software; if not, write to the Free Software Foundation, Inc., 51 Franklin
 Street, Fifth Floor, Boston, MA 02110-1301, USA 



1. Purpose of this software
----------------------------

The software here provided performs robust Stereo Odometry 
by using a novel technique called Coarse L1-Averaging.
Here we present a demo that calculate the trajectory of a moving
vehicle along a city. 

The user needs to provide feature associations for the different views, since
our software is just intended to estimate the trajectory. However, for the sake
of simplicity we have included the matcheds of a partial sequence from the KITTI
benchmark http://www.cvlibs.net/datasets/kitti/

For any doubt, comment, suggestion or bug just send me an e-mail to gros@cvc.uab.es

2. Software and execution
--------------------------------------
This software requires Matlab to run. Call the L1AVGVO function to run it.
Remember to set up all the parameters as explained below.

 This code performs Robust Stereo Visual Odometry from a set of
 matches between two views provided in Data and a calibration structure "calib".

 	Input:
		Data: 		Nx8 structure containing point matches (u, v) between two stereo views
				i.e., (left, right)_tk+1 (left, right)_tk

		calib:  	A structure containing calib.K a 3x3 matrix of intrinsic parameters
				and calib.B a scalar that specifies the baseline in meters

		options:	A structure with all the parameters required to run the method:
			* [options.nModels] Number of putative models generated from the input data (e.g., 100 - 1000)
			* [options.maxIters] Maximum amount of iterations for the iterative optimization of L1 averaging (e.g., 100- 1000)
			* [options.threshold] Projection threshold for an optional outlier refinement step (e.g., 0.5 - 3.5)	
			* [options.refinement] In case you want to activate the optional outlier refinement (true, false)
			* [options.subModels] Number of models selected for the averaging (the top K models) (e.g., 500)
			* [options.epsilon] Used in the stop condition of the L1 averaging (e.g., 1e-6)


	Output:
		sol:	A relative pose [R, T] in SE(3) as a 3x4 matrix


 This is a Matlab version of a C++ code available at https://github.com/germanRos/l1avgvo
 If you find any problem I would really appreciate if you inform me by e-mail.
 Of course, the C++ implementation is much faster, so in case you want to compare times, please
 always refer to the C++ version.


function [sol] = L1AVGVO(Data, calib, options, ~)#



3. Citing this software
------------------------

If you use this software for any of your work, please cite it as follows:

@inproceedings{ros13,
 author = {Ros, G. and Guerrero, J. and Sappa, A. D. and Ponsa, D. and L\'{o}pez, A. M.},
 title = {Fast and Robust l1-averaging-based Pose Estimation for Driving Scenarios},
 booktitle = {Proceedings of the British Machine Vision Conference},
 year = {2013},
 address = {Bristol, UK},
}#
