function fn =unscaledata(x,dmin,dmax)

fn = dmin + (dmax-dmin)*(x+1)/2;

end