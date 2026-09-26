import uproot
import awkward as ak
import numpy as np

f = uproot.open("/eos/user/a/arnaik/Skimmer_HtoAAto4g/test_input.root")
t = f["Events"]

photons = t.arrays(["Photon_pt", "Photon_eta", "Photon_pixelSeed"], library="ak")

pt  = photons["Photon_pt"]
eta = photons["Photon_eta"]
seed = photons["Photon_pixelSeed"]

n_all = len(pt)
print("total events:", n_all)

# Cut 1: >= 4 photons
nPho = ak.num(pt, axis=1)
mask = nPho >= 4
print("after >=4 photons:", int(ak.sum(mask)))

# leading 4, padded so <4-photon events don't crash (they're already False above)
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

print("FINAL survivors (no HLT):", int(ak.sum(mask)))
