import pickle as pkl
import ROOT

# only works for RobustScaler
scaler = pkl.load(open("../models/scaler_EB.pkl", "rb"))

# # extract the attributes
# center_values = scaler.center_
# scale_values = scaler.scale_

# # save the 
# data = {"center": center_values, "scale": scale_values}
# df = ROOT.RDF.MakeNumpyDataFrame(data)
# df.Snapshot("RobustScaler", "test.root")

a = [[i+1 for i in range(27)]]
print(scaler.transform(a))

ROOT.gInterpreter.ProcessLine(""" #include "../interface/MyRobustScaler.h" """)
scaler_test = ROOT.MyRobustScaler()
scaler_test.Load("test.root")
b = ROOT.std.vector("float")(i+1 for i in range(27))
print(scaler_test.Transform(b))


