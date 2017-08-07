function [yo] = UniPKLPE2(Y, A, po, va, di, vid, sid, y0, type)

%	This Code is based on the "Sect. 2.2	Least Power Estimation" and "Sect. 2.3 Spatially Correlated LPE" in the article 
%	"New Spatial Upscaling Methods for Multi-point Measurements: From Normal to p-Normal".
%	https://doi.org/10.1016/j.cageo.2017.08.001
%	Author: Feng Liu (unilf@126.com)
%	Date: 20170807
%
%   Y: Observation vector. 
%   A: Parameter Matrix with fields of lambda weight and calculated value like (0.038, 20.01).
%   Pow: p vector. 
%   va: variance vector. 
%   di, vid, sid: data flag (maybe they are not useful.)
%   y0: initial estimator (from least square estiamtion)
%
%   Computation of Lp Estimator with Model Y = AX + v
%
%   type: 0. LPE
%   type: 1. (unavailable) The method of IRLS(Iterative Reweighted Least Squares) .<1. Kuruoglu, 1998, Least Lp-norm impulsive noise  cancellation with polynomial filters. 2. Chen, 1994, Parameter estimation of linear systems with input-output noisy data A generalized lp norm approach.>
%   type: 2. (unavailable) The method of RSD(residual steepest descent).<1. Udahemuka, lp deconvolution for water pipe channel identification.>
%   type: 3. (unavailable) The method of Newton-Raphson. <Chen, 1994, Parameter estimation of linear systems with input-output noisy data A generalized lp norm approach.>
%   type: 4. GSLPE The direct method of Iterative gradient. <Nyquist, 1980, Recent studies on Lp-norm estimation.>
%   type: 5. Inverse distance weighting(IDW)


Pow = zeros(size(vid, 1), 1);
Var = zeros(size(vid, 1), 1);
Dis = zeros(size(vid, 1), 1);
nlist = [];
for i = 1 : size(vid, 1)
	lsid = vid(i);
    if lsid < 60
        nlist(size(nlist, 2) + 1) = i;
    end
	for j = 1 : size(sid, 1);
        if lsid == sid(j)
            Pow(i) = po(j);
            Var(i) = va(j);
            Dis(i) = di(j);
            break;
        end
	end
end

e = 0.01;
ve = 0.01;
Lambda = sqrt(gamma(3 ./ Pow) ./ gamma(1 ./ Pow));
SGM = ((Lambda ./ Var).^Pow);
Ws = (Pow) .* SGM;
icount = 100;

switch type
    case 0
        yo = 0;
        
        B = A(:, 1);
        C = A(:, 2);

        yls = y0;
        V = abs(Y - B * y0 - C);
        count = 0;
        while (abs(yo - y0) > ve && count <= icount)
            for i = 1 : size(V)
                if V(i) < e
                    V(i) = e;
                end
            end
            Wd = Ws .* (V.^(Pow - 2));
            W = diag(Wd);

            y0 = yo;
%             M = (B' * W * B) \ (B' * W);            yo = M * Y;
            yo =   (B' * W * B) \ (B' * W * Y);
            V = abs(Y - B * yo - C);
            count = count + 1;
        end
        if count > icount || isnan(yo) || (abs(yo - yls) > yls)
            yo = yls;
        end
    case 4
        B = A(:, 1);
        C = A(:, 2);
        
        coe = 1 ./ B;
        C = coe .* C;
        B = coe .* B;
        Yc = coe .* Y;

       yls = y0;
       V = abs(Yc - B * y0 - C);

       yo = 0;
       count = 0;
        while (abs(y0 - yo) > ve && count <= icount)
            for i = 1 : size(V)
                if V(i) < e
                    V(i) = e;
                end
            end
            Wd = Ws .* (V.^(Pow - 2));
            W = diag(Wd);

            y0 = yo;
            yo =   (B' * W * B) \ (B' * W * Y);
           V = abs(Yc - B * yo - C);

           count = count + 1;
        end
        
        if count >= icount || isnan(yo) || (abs(yo - yls) > yls)
            yo = yls;
        end
    case 5
        ds = zeros(size(Dis, 1), 1);
        s = 0.0;
        for i = 1 : size(Dis, 1)
            d = Dis(i);
            ds(i) = 1 / d^1;
            s = s + 1 / d^1;
        end
        ds = ds / sqrt(s);
        
        yo = 0;
        
        B = A(:, 1);
        C = A(:, 2);

        yls = y0;
        V = abs(Y - B * y0 - C);
        count = 0;
        while (abs(yo - y0) > ve && count <= icount)
            for i = 1 : size(V, 1)
                if V(i) < e
                    V(i) = e;
                end
            end
            Wd = Ws .* (V.^(Pow - 2));
            W = diag(ds .* Wd);

            y0 = yo;
%             M = (B' * W * B) \ (B' * W);            yo = M * Y;
            yo =   (B' * W * B) \ (B' * W * Y);
            V = abs(Y - B * yo - C);
            count = count + 1;
        end
        if count > icount || isnan(yo) || (abs(yo - yls) > yls)
            yo = yls;
        end
end




