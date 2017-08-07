function [Y, A] = UniPKModel(Data, h, to, radius, EPX, EPY)

%	This Code is based on the "Sect. 2.3 Spatially Correlated LPE" in the article 
%	"New Spatial Upscaling Methods for Multi-point Measurements: From Normal to p-Normal".
%	https://doi.org/10.1016/j.cageo.2017.08.001
%	Author: Feng Liu (unilf@126.com)
%	Date: 20170807
%
%   Data's fields are time, value, siteid, CoorX, CoorY, for example: [62217,12,120,618965.258602000,4301811.17807000;]
%   h is the step of semivariogram
%   to is the tolerance of semivariogram
%   radius is the search radius of a large-scale estimate
%   EPX, EPY: Geometrical center coordinate of a large-scale estimate
%
%   Computation of Model (Y=Ax+V) with Data, reverse Ordinary Kriging
%   Calculate the parameters Matrix A with Ordinary Kriging reversely.
%
%   We see the centre of the upscaling region as known point and make it as
%   a component of sample vector x. Whereas Y is the vector of observation
%   comprises all Zi.
%
%   step1: Calculate and fit the Semivariogram function with value, Coor and
%   radius 
%   step2: Calculate the weighted vector with Semivariogram and OK system
%   step3: Calculate the A with weighted vector and value

ColT = 1;
ColV = 2;
ColID = 3;
ColX = 4;
ColY = 5;

%step1
[xtemp, ytemp, ae] = Semivariogram(Data, h, to, radius);

%step2
%Build the matrix, calculate all semivariogram at first
dn = size(Data, 1);
UD = [];
for i = 1 : dn
    UD(i, :) = Data(i, :);
end
UD(dn + 1, :) = [Data(dn, ColT), 0, 0, EPX, EPY];

GMA = [];
for i = 1 : dn + 1
    x = [];
    for j = 1 : dn + 1
        d = sqrt((UD(i, ColX) - UD(j, ColX))^2 + (UD(i, ColY) - UD(j, ColY))^2);
        x(j) = d;
    end
    %g = Gaussian(ae, x);    %Exponential Spherical Gaussian
    g = Exponential(ae, x);
    GMA(i, :) = g;
end

A = [];
for i = 1 : dn
    % Build the OK system for each sample
    GM = GMA;
    GM(i, :) = [];
    GM0 = GM(:, i);
    GM0 = [GM0' 1]';
    GM(:, i) = [];
    b = zeros(dn, 1) - 1;
    GM = [GM b];
    b = zeros(1, dn) + 1;
    b = [b 0];
    GM = [GM' b'];
    GM = GM';
    
    if det(GM) == 0
        %break;
    end
    %Calculate parameter vector of lambda  GM * LD = GM0
    LD = GM \ GM0;
    
    % Calculate A
    VL = Data(:, ColV);
    VL(i, :) = [];
    sum = 0.0;
    for j = 1 : dn - 1
        sum = sum + VL(j) * LD(j);
    end
    A(i, :) = [LD(dn) sum];
end

Y = Data(:, ColV);
if det(GM) == 0
    %Y = [];    
end
end





