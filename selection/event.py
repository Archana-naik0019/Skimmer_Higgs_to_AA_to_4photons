import numpy as np
import awkward as ak

GAP_BARREL_ETA = 1.4442
GAP_ENDCAP_ETA = 1.5660
MAX_ETA = 2.5

def event_mask(
    collections,
    trigger_mask=None,
    apply_bJet_tagger=None,
    cut_4photons=True,
    cut_eta=True,
    cut_pixel_seed=True,
    cut_pt=True,
):
    photons = collections["Photon"]
    nPho = ak.num(photons, axis=1)

    mask = trigger_mask if trigger_mask is not None else ak.ones_like(nPho, dtype=bool)
    cut_masks = {}

    if cut_4photons:
        mask = mask & (nPho >= 4)
    cut_masks["4photons"] = mask

    photons_padded = ak.pad_none(photons, 4, axis=1)
    pho4 = photons_padded[:, :4]

    if cut_eta:
        abs_eta = np.abs(pho4.eta)
        eta_pass_each = (abs_eta < GAP_BARREL_ETA) | (
            (abs_eta > GAP_ENDCAP_ETA) & (abs_eta < MAX_ETA)
        )
        mask = mask & ak.fill_none(ak.all(eta_pass_each, axis=1), False)
    cut_masks["eta"] = mask

    if cut_pixel_seed:
        pixel_pass_each = pho4.pixelSeed == False
        mask = mask & ak.fill_none(ak.all(pixel_pass_each, axis=1), False)
    cut_masks["pixelSeed"] = mask

    if cut_pt:
        pt_pass_each = pho4.pt > 12.0
        mask = mask & ak.fill_none(ak.all(pt_pass_each, axis=1), False)
    cut_masks["pt_cuts"] = mask

    return mask, cut_masks