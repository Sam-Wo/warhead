"""Curated clinical-validation + patient-toxicity reference for the compounds
that recur at the top of the public screens.

Purpose (WARHEAD is an ADC-payload program): connect screen potency to (a) whether
the agent is clinically validated, (b) its dose-limiting / characteristic patient
toxicity and the DLT organ (feeds G6), and (c) whether the chemotype is an
established or candidate ADC payload.

Sources: FDA labels / standard oncology references / ADC literature. This is a
hand-curated table, not derived from the screens; treat toxicity entries as the
characteristic clinical pattern, not a per-patient prediction.
"""
from __future__ import annotations

import pandas as pd

# (compound, target/class, clinical status, approved example,
#  characteristic patient toxicity (DLT), DLT organ (G6 compartment), ADC-payload status)
_ROWS = [
    ("Romidepsin", "HDAC", "Approved", "CTCL / PTCL",
     "cytopenias, nausea, fatigue, QT prolongation", "bone marrow / cardiac", "not a payload"),
    ("Bortezomib", "Proteasome", "Approved", "multiple myeloma / MCL",
     "peripheral neuropathy, thrombocytopenia", "peripheral nerve / marrow", "not a payload"),
    ("Dactinomycin", "RNA Pol I (transcription)", "Approved", "Wilms, rhabdomyosarcoma",
     "myelosuppression, mucositis, hepatotoxicity (VOD)", "bone marrow / GI / liver", "not a payload"),
    ("Docetaxel", "Microtubule (stabiliser)", "Approved", "breast, NSCLC, prostate",
     "neutropenia, neuropathy, fluid retention", "bone marrow / peripheral nerve", "not a payload"),
    ("Paclitaxel", "Microtubule (stabiliser)", "Approved", "breast, ovarian, NSCLC",
     "neutropenia, peripheral neuropathy, hypersensitivity", "bone marrow / peripheral nerve", "not a payload"),
    ("Cabazitaxel", "Microtubule (stabiliser)", "Approved", "mCRPC (post-docetaxel)",
     "neutropenia (severe), diarrhoea", "bone marrow / GI", "not a payload"),
    ("Vinblastine", "Microtubule (destabiliser)", "Approved", "Hodgkin, testicular",
     "myelosuppression, neuropathy", "bone marrow / peripheral nerve", "not a payload"),
    ("Vinorelbine", "Microtubule (destabiliser)", "Approved", "NSCLC, breast",
     "neutropenia, neuropathy, phlebitis", "bone marrow / peripheral nerve", "not a payload"),
    ("Eribulin", "Microtubule (destabiliser)", "Approved", "breast, liposarcoma",
     "neutropenia, peripheral neuropathy", "bone marrow / peripheral nerve", "candidate ADC payload (MORAb-202)"),
    ("Dolastatin-10", "Microtubule (auristatin class)", "Phase 2 (parent)", "— (tox-limited)",
     "neutropenia, peripheral neuropathy", "bone marrow / peripheral nerve",
     "PARENT of auristatins MMAE/MMAF (payloads in Adcetris, Padcev, Polivy...)"),
    ("Monomethyl-auristatin-E (MMAE)", "Microtubule (auristatin)", "Approved (as ADC)", "brentuximab/enfortumab vedotin",
     "as ADC: neutropenia, peripheral neuropathy (MMAE), ocular (MMAF)", "peripheral nerve / marrow / cornea",
     "established ADC payload"),
    ("Maytansine (DM1/DM4)", "Microtubule (maytansinoid)", "Approved (as ADC)", "trastuzumab emtansine (T-DM1)",
     "as ADC: thrombocytopenia, hepatotoxicity, ocular (DM4)", "marrow / liver / cornea",
     "established ADC payload"),
    ("SN-38", "Topoisomerase I", "Approved (as irinotecan / ADC)", "irinotecan; sacituzumab govitecan",
     "diarrhoea (severe), neutropenia", "GI / bone marrow", "established ADC payload (Trodelvy)"),
    ("Exatecan / DXd", "Topoisomerase I", "Approved (as ADC)", "trastuzumab deruxtecan (Enhertu)",
     "as ADC: interstitial lung disease (ILD), neutropenia, nausea", "lung (ILD) / marrow / GI",
     "established ADC payload (deruxtecan)"),
    ("Camptothecin / Topotecan", "Topoisomerase I", "Approved (topotecan)", "SCLC, ovarian, cervical",
     "neutropenia, thrombocytopenia", "bone marrow", "TOP1 payload chemotype"),
    ("Nemorubicin (PNU-159682)", "Topoisomerase II / DNA", "Preclinical (metabolite)", "— (ADC payload)",
     "as ADC: haematologic (potent anthracycline)", "bone marrow", "candidate ADC payload (anthracycline)"),
    ("Epirubicin / anthracyclines", "Topoisomerase II / DNA", "Approved", "breast, many",
     "cardiotoxicity (cumulative), myelosuppression", "cardiac / bone marrow", "anthracycline payload chemotype"),
    ("Gemcitabine", "Ribonucleotide reductase / DNA", "Approved", "pancreatic, NSCLC, bladder",
     "myelosuppression, flu-like, transaminitis", "bone marrow / liver", "not a payload"),
    ("Bleomycin", "DNA strand breaks", "Approved", "Hodgkin, germ-cell",
     "pulmonary fibrosis (cumulative), skin", "lung", "not a payload"),
    ("Epothilone-b / Ixabepilone", "Microtubule (stabiliser)", "Approved (ixabepilone)", "breast",
     "peripheral neuropathy, neutropenia", "peripheral nerve / marrow", "candidate payload chemotype"),
    ("Triptolide / Minnelide", "XPB/ERCC3 (transcription)", "Phase 1/2 (prodrug)", "— (investigational)",
     "GI toxicity, myelosuppression (narrow window)", "GI / bone marrow", "explored as payload"),
    ("Sepantronium (YM155)", "Survivin (BIRC5) suppressant", "Phase 2 (discontinued)", "—",
     "well tolerated; febrile neutropenia in some", "bone marrow", "not a payload"),
    ("Daporinad (FK866/APO866)", "NAMPT", "Phase 1/2 (discontinued)", "—",
     "thrombocytopenia, lymphopenia", "bone marrow", "candidate payload chemotype (NAMPTi)"),
    ("Trametinib", "MEK1/2", "Approved", "BRAF-mut melanoma, NSCLC (combo)",
     "rash, diarrhoea, cardiomyopathy, ocular (RVO)", "skin / GI / cardiac / retina", "not a payload"),
    ("Binimetinib / Pimasertib", "MEK1/2", "Approved (binimetinib)", "BRAF-mut melanoma (combo)",
     "rash, CK elevation, retinopathy", "skin / muscle / retina", "not a payload"),
    ("Dinaciclib", "CDK1/2/5/9", "Phase 3 (discontinued)", "—",
     "myelosuppression, GI", "bone marrow / GI", "not a payload"),
    ("Luminespib (AUY922)", "HSP90", "Phase 2", "—",
     "diarrhoea, night blindness (ocular), fatigue", "GI / retina", "not a payload"),
    ("Pralatrexate", "DHFR / antifolate", "Approved", "PTCL",
     "mucositis, thrombocytopenia", "GI / bone marrow", "not a payload"),
    ("Staurosporine", "Pan-kinase (tool)", "Preclinical (tool)", "—",
     "n/a (not clinically used)", "n/a", "not a payload (tool compound)"),
]

COLUMNS = ["compound", "target_class", "clinical_status", "approved_example",
           "patient_toxicity_DLT", "dlt_organ", "adc_payload_status"]


def clinical_tox_table() -> pd.DataFrame:
    return pd.DataFrame(_ROWS, columns=COLUMNS)
