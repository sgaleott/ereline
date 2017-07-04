import pandas as pd
import sys

filename = sys.argv[1]
chtag = filename.split("_")[-1].split(".")[0]

g = pd.read_hdf(filename.format(chtag),"data" )
gsm = smooth(g, chtag)
gsm.to_hdf(filename.replace(chtag, "smoothed_" + chtag), "data" )
