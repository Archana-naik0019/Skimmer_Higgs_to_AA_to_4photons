import numpy as np
import awkward as ak

# Detector geometry constants
GAP_BARREL_ETA = 1.4442
GAP_ENDCAP_ETA = 1.5660
MAX_ETA = 2.5

def photon_mask(photons, apply_pixelSeed=False):
    """
    Per-photon object-level mask. Applied inline in reducer.run()
    as obj[photon_mask(obj, apply_pixelSeed=self.apply_pixelSeed)].
    """
    abs_eta = np.abs(photons.eta)
    eta_pass = (abs_eta < GAP_BARREL_ETA) | (
        (abs_eta > GAP_ENDCAP_ETA) & (abs_eta < MAX_ETA)
    )

    pt_pass = photons.pt > 12.0

    mask = eta_pass & pt_pass

    if apply_pixelSeed:
        mask = mask & (photons.pixelSeed == False)

    return mask