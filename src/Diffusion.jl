module Diffusion
using StaticArrays 

export absorb_boundry
export reflect_boundary 

function absorb_boundry(pos,r)
    if(((pos[1]^2)+(pos[2]^2))^(1/2) < r)
    return pos
    else
    return nothing
    end
end

function reflect_boundary(oldpos,pos,r)
    if(((pos[1]^2)+(pos[2]^2))^(1/2) < r)
    return pos
    else
        
    end
end

end
