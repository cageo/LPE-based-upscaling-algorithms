function g = Exponential(x,hdata)%x=[range: a sill: c nugget: c0]
    g = x(3) + x(2) - x(2) * exp(-3 * hdata / x(1));