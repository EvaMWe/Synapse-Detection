%calculates the second order intensity moment

function parameter = m2Calc(parameter)
ij2 = parameter.ij2;
area = parameter.area;
m0 = parameter.m0;

m2 = (sum(sum(ij2.*area)))/m0;
parameter.m2 = m2;
end