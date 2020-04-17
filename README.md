# Orthogonal-Dist-to-Ellipse
MATLAB codes to compute the orthogonal distance between a point and an ellipse, and the orthogonal contacting point on the ellipse.

Three algorithms are implemented whose source code is deposited here: the algorithm proposed in (Ahn 2001) "Least-squares orthogonal distances fitting of circle, sphere, ellipse, hyperbola, and parabola," Pattern Recognition; the exact algorithm and the convergent iterative algorithm proposed in April, 2020 by Siyu Guo et al. 
The M-file name for the algorithms are ellipse_orthogonal_dist_arw01.m, ellipse_orthogonal_dist_exact.m and ellipse_orthogonal_dist_ci.m, respectively. 
Two C++ source files are also provided, namely, arw01_mex.cpp and ci_mex.cpp. Use mex to compile them into the MEX modules used in the M-files.
exprm_1.m, exprm_2.m and exprm_3.m are three demonstrations without input or output. draw_ellipse.m is used in the demonstrations to draw an ellipse in a MATLAB figure.
