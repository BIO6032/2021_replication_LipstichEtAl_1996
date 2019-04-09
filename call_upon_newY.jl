include("Fig1_strain_parameters.jl")

function thousand_para()
    ux = 0.2;
    for k in 1:1000
        array = zeros((1000,4))
        b = new_Y()
        array[k,:] = [b.uy, b.by, b.ey, b.βy]
        return array
    end
end

thousand_para()
