# Convert netCDF file to ancil using ants

import iris, sys
import ants.io.save as save

cube = iris.load_cube(sys.argv[1])
save.ancil(cube, sys.argv[2])

