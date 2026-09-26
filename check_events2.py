import uproot
import awkward as ak
import numpy as np

f = uproot.open("/eos/user/a/arnaik/Skimmer_HtoAAto4g/test_input.root")
t = f["Events"]

data = t.arrays(
    ["Photon_pt", "Photon_eta", "Photon_pixelSeed",
     "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"],
    library="ak"
)

pt   = data["Photon_pt"]
eta  = data["Photon_eta"]
seed = data["Photon_pixelSeed"]
hlt  = data["HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"]

n_all = len(pt)
print("total events:", n_all)

# HLT first, matching your reducer's order
mask = hlt
print("after HLT:", int(ak.sum(mask)))

# Cut 1: >= 4 photons
nPho = ak.num(pt, axis=1)
mask = mask & (nPho >= 4)
print("after >=4 photons:", int(ak.sum(mask)))

pt4  = ak.pad_none(pt, 4, axis=1)[:, :4]
eta4 = ak.pad_none(eta, 4, axis=1)[:, :4]
seed4 = ak.pad_none(seed, 4, axis=1)[:, :4]

GAP_BARREL_ETA = 1.4442
GAP_ENDCAP_ETA = 1.5660
MAX_ETA = 2.5

abs_eta = np.abs(eta4)
eta_pass_each = (abs_eta < GAP_BARREL_ETA) | ((abs_eta > GAP_ENDCAP_ETA) & (abs_eta < MAX_ETA))
eta_pass_event = ak.fill_none(ak.all(eta_pass_each, axis=1), False)
mask = mask & eta_pass_event
print("after eta cut:", int(ak.sum(mask)))

pixel_pass_each = seed4 == False
pixel_pass_event = ak.fill_none(ak.all(pixel_pass_each, axis=1), False)
mask = mask & pixel_pass_event
print("after pixelSeed cut:", int(ak.sum(mask)))

pt_pass_each = pt4 > 12.0
pt_pass_event = ak.fill_none(ak.all(pt_pass_each, axis=1), False)
mask = mask & pt_pass_event
print("after pt cut:", int(ak.sum(mask)))

print("FINAL survivors (with HLT):", int(ak.sum(mask)))
