function [x, y, ae] = Semivariogram(v, h, to, radius)

%	Calculate the Semivariogram
%
%   v is the data, its fields are time, value, siteid, CoorX, CoorY, for example: [62217,12,120,618965.258602000,4301811.17807000;]
%   h is the step of semivariogram
%   to is the tolerance of semivariogram
%   radius is the search radius of a large-scale estimate
%
%

Colid = 1;
Collon = 2;
Collat = 3;
ColX = 4;
ColY = 5;

ColTime = 1;
ColVal = 2;

%range
x = [h: h: radius*2];
%x = [h: h/2: radius radius+h: h: radius*2];
y = zeros(1,size(x,2));
ycount = zeros(1,size(x,2));

sum = 0.0;
for i = 1 : size(v, 1) - 1
    Xi = v(i, ColX);
    Yi = v(i, ColY);
    Vi = v(i, ColVal);
    
    for j = i + 1 : size(v, 1)
        Xj = v(j, ColX);
        Yj = v(j, ColY);
        Vj = v(j, ColVal);
        
        d = sqrt((Xi - Xj)^2 + (Yi - Yj)^2);
        re = roundn(mod(d, h), -4);
        if re >= (h-to) || re <= to
        %if re >= (h-to) || re <= (h+to)
            idx = round(d / h);
            if idx <= numel(y) && idx > 0
                y(idx) = y(idx) + (Vi - Vj)^2;
                ycount(idx) = ycount(idx) + 1;
            end
        end
    end
end

y = y ./ ycount;

ys = size(y, 2);
for i = ys : -1 : 1
    if isnan(y(i)) || (i > ys / 2 && mod(i, 4) ~= 1)
        y(i) = [];
        x(i) = [];
    end
end


x0 = [1 1 1] ;
% Only the Exponential model was considered
ae = lsqcurvefit('Exponential', x0, x, y);
fe = Exponential(ae, x);

end

