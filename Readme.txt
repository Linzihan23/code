***********************************************************************************************************
***********************************************************************************************************

Matlab demo code for "Digital holographic microscopy phase noise reduction based on an over-complete
chunked discrete cosine transform sparse dictionary" 

by£ºZihan Lin (linzihan21@stu.xjtu.edu.cn)

If you use/adapt our code in your work (either as a stand-alone tool or as a component of any algorithm),
you need to appropriately cite our paper.

This code is for academic purpose only. Not for commercial/industrial activities.

Note"
The runtimes reported in the paper are from a lab computer implementation. This Matlab version is intended 
to facilitate understanding of the algorithm. This code has not been optimized and its speed is not representative. 
The results may differ slightly from those in the paper due to cross-platform transfer.

***********************************************************************************************************
***********************************************************************************************************

Usage:

demo_main.m:Use the demo_main.m program to obtain a digital holographic reconstruction of the phase data map.

demo_DCT.m:Noise reduction experiments are performed using the demo_DCT.m program to simulate the method proposed in this paper.

***********************************************************************************************************
***********************************************************************************************************

Evaluation index:

ENL.m is the equivalent visual number evaluation index procedure.

SSI.m is the scatter suppression index evaluation index procedure.

SSIM.m is the structural similarity index evaluation index procedure.

NC.m is the normalized correlation coefficient evaluation index program.

FeatureSIM.m is the feature similarity index subroutine evaluation index program.

***********************************************************************************************************
***********************************************************************************************************

Simulation and test data:

guang.bmp: used to simulate experimental data graphs.

Interdigital electrode.tif: for testing and comparing experimental data graphs.