# Binary_FFPs_microlensing


There are several codes that generate  microlensing light curves due to binary free-floating exoplanets, with the names 1_Cteta.py, 2_Cu0.py, 
3_Crho.py, .... In these light curves, one parameter changes 23 times and each time new light curve is made.  Then these lightcurves using 
the following command line convert to *.mp4  file.  The command line is:  

cat LIc_0.jpg  LIc_1.jpg LIc_2.jpg  LIc_3.jpg LIc_4.jpg LIc_5.jpg LIc_6.jpg LIc_7.jpg LIc_8.jpg LIc_9.jpg LIc_10.jpg LIc_11.jpg  LIc_12.jpg  
LIc_13.jpg LIc_14.jpg LIc_15.jpg LIc_16.jpg  LIc_17.jpg LIc_18.jpg LIc_19.jpg LIc_20.jpg LIc_21.jpg LIc_22.jpg LIc_23.jpg| ffmpeg -f image2pipe
-r 1 -vcodec mjpeg -i - -vcodec  libx264 central_theta.mp4


All *.mp4 files were made according to these codes can be found in  the following link:  
https://iutbox.iut.ac.ir/index.php/s/NYnnowMoLbFnZeD


The C++ codes (*.cpp) do Monte Carlo simulations from microlensing events due to binary free-floating exoplanets. These codes use 
VBBinaryLensing package to calculate magnification factor for binary lensing from extended source stars.  

Codes *_stat.py  calculates the statistical properties from microlensing events due to binary free-floating planets.  

BFFs_case.cpp makes a bis ensamble of all possible microlensing events due to each known moon-planet system orbiting the Sun. 

"dist_BFFPs.py" plots the (q, and d) distributions of the known moon-planet systems.  
