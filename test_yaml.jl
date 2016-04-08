using YAML
using Dolo

data = YAML.load_file("/home/pablo/Programming/econforge/dolo/examples/models/rbc.yaml")
# data = YAML.load_file("rbc.yaml")
data["equations"]

recipe = RECIPES[:dtcscc]

methods(SymbolicModel)

model = SymbolicModel(data, :dtcscc, "rbc.yaml")


model.calibration[:z]
calib = ModelCalibration(model)
calib.flat

calib.grouped[:states]


nmodel = DTCSCCModel(model)

s,x,y,e,p = nmodel.calibration[:states,:controls,:auxiliaries,:shocks,:parameters]
f = nmodel.functions.arbitrage
methods(evaluate)
f = nmodel.functions.arbitrage
g = nmodel.functions.transition
v = nmodel.functions.value
a = nmodel.functions.auxiliary
#
# f(s,x,y,e,s,x,y,p)
model.equations

funs = nmodel.functions


Y = evaluate(a,s,x,p)
S = evaluate(g,s,x,e,p)

res = evaluate(f,s,x,e,s,x,p)



N = 100


S = repmat(s',N)
X = repmat(x',N)
E = repmat(e',N)
P = repmat(p',N)

evaluate(f,S,X,E,S,X,P)
