using InstrumentOperator
using NCDatasets
using OrderedCollections
using YAML

# Load L1 file (replace with your example!):
L1File = "/net/fluo/data1/group/oco2/L1bSc/oco2_L1bScND_26780a_190715_B10003r_200429212407.h5"
ocoData = Dataset(L1File);

# Load corresponding met file:
metFile = "/net/fluo/data1/group/oco2/L2Met/oco2_L2MetND_26780a_190715_B10003r_200429212406.h5"
metData = Dataset(metFile);

# Load dictionary:
dictOCO2 = YAML.load_file("json/oco2.yaml"; dicttype=OrderedDict{String,Any});

# Load L1 file (could just use filenames here as well)
oco = InstrumentOperator.load_L1(dictOCO2,ocoData, metData);

# Pick some bands as tuple (or just one)
bands = (1,2);
# Indices within that band:
indices = (1:1016,20:1000);
# Geo Index (footprint,sounding):
GeoInd = [5,3000];

# Get data for that sounding:
oco_sounding = InstrumentOperator.getMeasurement(oco, bands, indices, GeoInd)