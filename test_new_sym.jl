
import Dolo


url = "examples/models/rbc_dtcc_mc.yaml"
typeof(url)

model = Dolo.Model(url)
data = model.data

@time Dolo.get_name(model)
@time Dolo.get_symbols(model)
@time Dolo.get_definitions(model)
@time Dolo.get_equations(model)
@time Dolo.get_infos(model)
@time Dolo.get_options(model)
@time Dolo.get_definitions(model)

# apparently the time is dominated by the solution of the calibration
# the result should be cached

@time calib = Dolo.get_calibration(model)



syms


Dolo.ModelCalibration(calibs, syms)





@time domain = NewDolo.get_domain(model)
@time exo = NewDolo.get_exogenous(model)
@time grid = NewDolo.get_grid(model)


calib = NewDolo.get_calibration(model)
syms = NewDolo.get_symbols(model)
