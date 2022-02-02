# Computing some distributions for discrete values
using Distributions

function make_distributions(c_min,c_max,number,N)

values = linspace(c_min,c_max,number)

    return sample(values,N)

end
