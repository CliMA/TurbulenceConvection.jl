
function center_field(grid::Grid)
    return pyzeros(grid.nzg)
end

# TODO: change to use unequal length to center_field
function face_field(grid::Grid)
    return pyzeros(grid.nzg)
end
