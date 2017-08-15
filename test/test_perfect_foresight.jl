import Dolo

path = Dolo.pkg_path

@testset "testing perfect_foresight" begin

    fn = joinpath(path, "examples", "models", "rbc_dtcc_ar1.yaml")

    model = Dolo.yaml_import(fn)

    # we must define the series for exogenous shocks
    n_e = length(model.symbols[:exogenous])

    T_e = 5
    exo = zeros(T_e, n_e)
    exo[1,:] = [0.00]  # this is used to determine initial steady-state of the model
    exo[2,:] = [0.01]
    exo[3,:] = [0.02]
    exo[4,:] = [0.03]
    exo[5,:] = [0.04]  # this is used to determine final steady-state

    @time df = Dolo.perfect_foresight(model, exo, T=20)

    @test true

end


# Second method: Define a dictionary with shocks
# Compare the results with the previous method

path = Dolo.pkg_path
fn = joinpath(path, "examples", "models", "rbc_dtcc_ar1.yaml");

model = Dolo.yaml_import(fn);

n_e = length(model.symbols[:exogenous])
T_e = 5
exo = zeros(T_e, n_e)
exo[1,:] = [0.00]  # this is used to determine initial steady-state of the model
exo[2,:] = [0.01]
exo[3,:] = [0.02]
exo[4,:] = [0.03]
exo[5,:] = [0.04]
@time df_1matrix = Dolo.perfect_foresight(model, exo, T=20)

exo = Dict(:z =>  [0, 0.01, 0.02, 0.03, 0.04])
@time df_1dict = Dolo.perfect_foresight(model, exo, T=20)

if df_1dict == df_1matrix
  println("df is the same")
end


# Two shocks
fn = joinpath(path, "examples", "models", "rbc_catastrophe.yaml");

model = Dolo.yaml_import(fn)

n_e = length(model.symbols[:exogenous])
T_e = 5
exo = zeros(T_e, n_e)
exo[1,:] = [0.00, 0.00]  # this is used to determine initial steady-state of the model
exo[2,:] = [0.01, 0.00]
exo[3,:] = [0.02, 0.00]
exo[4,:] = [0.03, 0.00]
exo[5,:] = [0.04, 0.00]
@time  df_2matrix = Dolo.perfect_foresight(model, exo, T=20)


exo = Dict(:z =>  [0.00, 0.01, 0.02, 0.03, 0.04], :xi =>  [0.00, 0.00, 0.00, 0.00, 0.00])
@time  df_2dict_a = Dolo.perfect_foresight(model, exo, T=20)

if df_2matrix == df_2dict_a
  println("df is the same")
end

# or

exo = Dict(:z =>  [0.00, 0.01, 0.02, 0.03, 0.04])
@time df_2dict_b = Dolo.perfect_foresight(model, exo, T=20)

if df_2matrix == df_2dict_b
  println("df is the same")
end


#or order changed
exo = Dict( :xi =>  [0.00, 0.00, 0.00, 0.00, 0.00], :z =>  [0.00, 0.01, 0.02, 0.03, 0.04])
@time df_2dict_c = Dolo.perfect_foresight(model, exo, T=20)

if df_2matrix == df_2dict_c
  println("df is the same")
end

#or not defined as Float64
exo = Dict( :xi =>  [0, 0, 0, 0, 0], :z =>  [0.00, 0.01, 0.02, 0.03, 0.04])

@time df_2dict_d = Dolo.perfect_foresight(model, exo, T=20)
if df_2matrix == df_2dict_d
  println("df is the same")
end
