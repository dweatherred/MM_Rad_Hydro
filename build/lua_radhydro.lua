
timing = {}
timing["dt"] = 1.0E-12 
timing["time_initial"] = 0
timing["time_final"] = 5.0E-9  --seconds

spatial_cells = {}
spatial_cells["xl_bound"] = 0 --cm
spatial_cells["xr_bound"] = 2 -- cm
spatial_cells["num_cells"] = 9000

-- Shock wave input
--[[
timing = {}
timing["dt"] = 1.0E-12 
timing["time_initial"] = 0
timing["time_final"] = 3.0E-9  --seconds

spatial_cells = {}
spatial_cells["xl_bound"] = 0 --cm
spatial_cells["xr_bound"] = 1 -- cm
spatial_cells["num_cells"] = 9000
]]--

-- Sod Shock Input
--[[
timing = {}
timing["dt"] = 9.95586E-06
timing["time_initial"] = 0
timing["time_final"] = 0.2 --seconds

spatial_cells = {}
spatial_cells["xl_bound"] = 0 --cm
spatial_cells["xr_bound"] = 1 --cm
spatial_cells["num_cells"] = 1000
]] --

-- Marshak Input

--[[
timing = {}
timing["dt"] = 1.0E-12
timing["time_initial"] = 0
timing["time_final"] = 5.0E-8 --seconds

spatial_cells = {}
spatial_cells["xl_bound"] = 0 --cm
spatial_cells["xr_bound"] = 2 --cm
spatial_cells["num_cells"] = 80
]] --