# LPE-based-upscaling-algorithms
This is the core code of LPE, GSLPE and IDWLPE in the article "New Spatial Upscaling Methods for Multi-point Measurements: From Normal to p-Normal". https://doi.org/10.1016/j.cageo.2017.08.001
Authors: Feng Liu (unilf@126.com), Xin Li
All the codes were written in MATLAB. Since the M-files are formed as MATLAB functions, they cannot run in your computers directly.

Example:
  % If all the parameters are defined in advance
	% GSLPE
    [Y, A] = UniPKModel(v, h, to, radius, Xc, Yc);
    xlpk = UniPKLPE2(Y, A, Pow, Var, Dis, VID, SID, xls, 4);
    
  % LPE
    A = zeros(size(Y, 1), 2);
    A(:, 1) = A(:, 1) + 1;
    xlp = UniPKLPE2(Y, A, Pow, Var, Dis, VID, SID, xls, 0);
    
  % IDWLPE
    xlpi = UniPKLPE2(Y, A, Pow, Var, Dis, VID, SID, xls, 5);
