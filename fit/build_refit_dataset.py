"""Build the arm-level refit data files from the checked trial extraction.

Canonical builder from 2026-10-25 on. It replaces
References_refit_trials/search/extraction/build_refit_dataset.py, which built the
10-22 to 10-25 files and is kept unchanged as their record.

Usage (from the repository root):
    python3 -I fit/build_refit_dataset.py OUT.csv [flags]
Extraction folder: env BED_EXTRACT_DIR, default ../References_refit_trials/search/extraction.

Main set (Igor Shuryak, 2026-10-22 and 2026-10-25): definitive radiotherapy alone of
non-nasopharyngeal HNSCC, randomized, the 2018 Table 1 trials re-extracted from their papers,
ARTSCAN and the three trials added on 2026-10-24 in sensitivity sets only.

Changes on 2026-10-25 (review of 2026-10-25, plan decisions 4-7):
- KROG HYPO arm: 65 T1 + 9 T2 patients (Moon 2014 Table 1: T2a 7 + T2b 2), not 8.
- RTOG 9003: rates and patient numbers from the same report (Fu 2000 Table 2 and Table 1,
  2-year KM local-regional control). Beitler 2014 updated rates are a sensitivity set
  (--rtog9003-beitler); its per-arm N is not reported, so it borrows the 2000 N.
- Mixed fraction sizes are kept as components: cell-kill BED = sum(d) + sum(d^2)/(alpha/beta),
  written as total_dose_gy and sum_d2_gy2.
- Overall time: the source time (delivered median or mean, else planned) stays the main input.
  T_nominal_days is the planned schedule on one calendar (first fraction on a Monday, last
  fraction day inclusive), for the nominal-clock sensitivity set (--nominal-clock).
- The estimator, the endpoint definition and the exact source percentage are carried into the file.
- The builder refuses rows whose patient number differs from the matched extraction row
  without a stated reason, and fraction or patient components that do not sum.
"""
import csv, json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTRACT = os.environ.get("BED_EXTRACT_DIR", os.path.join(ROOT, "..", "References_refit_trials", "search", "extraction"))

MON_FRI = (1, 1, 1, 1, 1, 0, 0)
MON_SAT = (1, 1, 1, 1, 1, 1, 0)
EVERY_DAY = (1, 1, 1, 1, 1, 1, 1)
BID_MON_FRI = (2, 2, 2, 2, 2, 0, 0)


def last_day(F, week=MON_FRI, gap=None):
    """Inclusive day of the last fraction, first fraction on Monday (day 1).
    week: fractions given on Mon..Sun. gap: (after this many fractions, rest this many calendar days)."""
    done, day, rested = 0, 0, 0
    while done < F:
        day += 1
        if gap and done >= gap[0] and rested < gap[1]:
            rested += 1
            continue
        done += week[(day - 1) % 7]
    return day


def nominal(spec):
    """spec: list of (weight, F, week, gap). Weighted mean of last_day; a non-integer F is
    interpolated between its floor and ceiling."""
    tot = w = 0.0
    for wt, F, week, gap in spec:
        lo, hi = math.floor(F), math.ceil(F)
        t = last_day(lo, week, gap) if lo == hi else \
            (hi - F) * last_day(lo, week, gap) + (F - lo) * last_day(hi, week, gap)
        tot += wt * t; w += wt
    return round(tot / w, 3)


def one(F, week=MON_FRI, gap=None):
    return [(1, F, week, gap)]


# unit, trial, pmid, arm match (substring of extracted arm label), d (mean per fraction), nfx, T (source),
# N, value (as printed), failure (True if value is a failure %), T rule, notes
A = [
 ("CAIR_2.0Gy", "CAIR", "17011446", "CAIR, 2.0 Gy", 2.0, 35, 36, 23, 73, False, "median delivered", "2.0 Gy period (1993-94)"),
 ("CAIR_2.0Gy", "CAIR", "17011446", "Control, 2.0 Gy", 2.0, 35, 51, 21, 30, False, "median delivered", "2.0 Gy period (1993-94)"),
 ("CAIR_1.8Gy", "CAIR", "17011446", "CAIR, 1.8 Gy", 1.8, 40, 40, 28, 78, False, "median delivered", "1.8 Gy period (1995-96)"),
 ("CAIR_1.8Gy", "CAIR", "17011446", "Control, 1.8 Gy", 1.8, 40, 54, 28, 35, False, "median delivered", "1.8 Gy period (1995-96)"),
 ("CAIR_vs_CB", "CAIR_vs_CB", "22836063", "CAIR", 1.8, 69.1 / 1.8, 39, 172, 63, False, "mean delivered", "fractions = mean delivered dose 69.1 Gy / 1.8 Gy"),
 ("CAIR_vs_CB", "CAIR_vs_CB", "22836063", "CB", 1.8, 70.1 / 1.8, 40, 171, 65, False, "mean delivered", "fractions = mean delivered dose 70.1 Gy / 1.8 Gy"),
 ("CHART", "CHART", "20394851", "CHART", 1.5, 36, 12, 550, 45, False, "planned", "5-y value read from Fig. 3a"),
 ("CHART", "CHART", "20394851", "Conventional", 2.0, 33, 45.5, 365, 50, False, "planned", "5-y value read from Fig. 3a"),
 ("DAHANCA67", "DAHANCA67", "14511925", "6 fractions", 2.0, 33.5, 39, 750, 70, False, "median delivered", "33-34 fractions"),
 ("DAHANCA67", "DAHANCA67", "14511925", "5 fractions", 2.0, 33.5, 46, 726, 60, False, "median delivered", "33-34 fractions"),
 ("EORTC22791", "EORTC22791", "1480768", "CF", 2.0, 35, 51, 159, 40, False, "median delivered", ""),
 ("EORTC22791", "EORTC22791", "1480768", "HF", 1.15, 70, 51, 166, 59, False, "median delivered", ""),
 ("EORTC22811", "EORTC22811", "7569030", "SDF", 2.0, 35, 49, 168, 27, False, "planned", "arm allowed 70 Gy/35/7 wk or 75 Gy/44/9 wk, split not reported; 70 Gy used; 5-y value read from Fig. 3"),
 ("EORTC22811", "EORTC22811", "7569030", "MFD (three", 1.6, 43.5, 45.5, 163, 23, False, "planned (6-7 wk stated)", "42 or 45 fractions; misonidazole arm excluded; 5-y value read from Fig. 3"),
 ("EORTC22851", "EORTC22851", "9288839", "CF", 2.0, 35, 54, 255, 46, False, "median delivered", ""),
 ("EORTC22851", "EORTC22851", "9288839", "AF", 1.6, 45, 33, 257, 59, False, "median delivered", ""),
 ("Ghoshal", "Ghoshal", "18343310", "Concomitant boost", 67.5 / 40, 40, 36, 143, 73.6, False, "median delivered", "45 Gy/25 + 22.5 Gy/15 boost"),
 ("Ghoshal", "Ghoshal", "18343310", "Conventional", 2.0, 33, 46, 142, 54.5, False, "median delivered", ""),
 ("Hliniak", "Hliniak", "11830307", "CF", 2.0, 33, 45, 193, 76, False, "planned", ""),
 ("Hliniak", "Hliniak", "11830307", "AF", 2.0, 33, 38, 188, 81, False, "planned", ""),
 ("IAEA_ACC", "IAEA_ACC", "20382075", "six fractions", 2.0, 34, 40, 452, 42, False, "median delivered", "33-35 fractions"),
 ("IAEA_ACC", "IAEA_ACC", "20382075", "five fractions", 2.0, 34, 47, 448, 30, False, "median delivered", "33-35 fractions"),
 ("Jackson", "Jackson", "9165135", "Once daily", 2.0, 33, 47.6, 41, 44.3, False, "mean delivered", ""),
 ("Jackson", "Jackson", "9165135", "Twice daily", 2.0, 33, 24.25, 41, 49.1, False, "mean delivered", ""),
 ("KROG0201", "KROG0201", "24161568", "CONV", 2.0, (74 * 33 + 8 * 35) / 82, 50, 82, 77.8, False, "median delivered", "fractions weighted by T1 74 / T2 8 (Table 1)"),
 ("KROG0201", "KROG0201", "24161568", "HYPO", 2.25, (65 * 28 + 9 * 30) / 74, 42, 74, 88.5, False, "median delivered", "fractions weighted by T1 65 / T2 9 (Table 1: T2a 7, T2b 2)"),
 ("PMH_HF", "PMH_HF", "32055717", "HF", 1.45, 40, 29, 169, 49, False, "median delivered", "competing-risk estimate"),
 ("PMH_HF", "PMH_HF", "32055717", "CF", 2.55, 20, 29, 162, 40, False, "median delivered", "competing-risk estimate"),
 ("RTOG8313_A", "RTOG8313", "2180866", "67.2", 1.2, 56, 39.2, 63, 25, False, "derived: 5 d/wk, 2 fx/d", "randomization scheme A only (Cox 1990)"),
 ("RTOG8313_A", "RTOG8313", "2180866", "72.0", 1.2, 60, 42.0, 58, 37, False, "derived: 5 d/wk, 2 fx/d", "randomization scheme A only (Cox 1990)"),
 ("RTOG8313_A", "RTOG8313", "2180866", "76.8", 1.2, 64, 44.8, 116, 42, False, "derived: 5 d/wk, 2 fx/d", "randomization scheme A only (Cox 1990)"),
 ("RTOG8809", "RTOG8809", "7790243", "AHFX-S", 1.6, 42, 43.5, 38, 45, True, "median delivered", "2-y local-regional failure 45% -> control 55%"),
 ("RTOG8809", "RTOG8809", "7790243", "AFX-C", 70.5 / 41, 41, 43, 32, 56, True, "median delivered", "2-y local-regional failure 56% -> control 44%"),
 ("RTOG9003", "RTOG9003", "10924966", "Standard", 2.0, 35, 50, 268, 46.0, False, "median delivered", "rate and N from Fu 2000 (Tables 2 and 1)"),
 ("RTOG9003", "RTOG9003", "10924966", "Hyperfractionation", 1.2, 68, 50, 263, 54.4, False, "median delivered", "rate and N from Fu 2000 (Tables 2 and 1)"),
 ("RTOG9003", "RTOG9003", "10924966", "with split", 1.6, 42, 43, 274, 47.5, False, "median delivered", "rate and N from Fu 2000 (Tables 2 and 1)"),
 ("RTOG9003", "RTOG9003", "10924966", "concomitant boost", 72 / 42, 42, 43, 268, 54.5, False, "median delivered", "rate and N from Fu 2000 (Tables 2 and 1)"),
 ("RTOG9512", "RTOG9512", "25035199", "SFX", 2.0, 35, 49, 119, 67, False, "planned", "T2N0 glottic"),
 ("RTOG9512", "RTOG9512", "25035199", "HFX", 1.2, 66, 45.5, 120, 73, False, "planned", "T2N0 glottic"),
 ("TROG9101", "TROG9101", "11439206", "CRT", 2.0, 35, 50, 171, 47, False, "median delivered", ""),
 ("TROG9101", "TROG9101", "11439206", "ART", 1.8, 33, 24, 172, 52, False, "median delivered", ""),
 ("Yamazaki_small", "Yamazaki", "16169681", "Arm A-1", 2.0, 30, 42, 31, 76, False, "planned", "tumor <= 2/3 of cord; T1 glottic"),
 ("Yamazaki_small", "Yamazaki", "16169681", "Arm B-1", 2.25, 25, 35, 31, 91, False, "planned", "tumor <= 2/3 of cord; T1 glottic"),
 ("Yamazaki_large", "Yamazaki", "16169681", "Arm A-2", 2.0, 33, 46.2, 57, 78, False, "planned", "tumor > 2/3 of cord; T1 glottic"),
 ("Yamazaki_large", "Yamazaki", "16169681", "Arm B-2", 2.25, 28, 39.2, 61, 92, False, "planned", "tumor > 2/3 of cord; T1 glottic"),
]
ARTSCAN = [
 ("ARTSCAN", "ARTSCAN", "26427805", "CF (conventional", 2.0, 34, 46, 367, 64.9, False, "planned", ""),
 ("ARTSCAN", "ARTSCAN", "26427805", "AF (accelerated", 68 / 43, 43, 31, 366, 65.5, False, "planned", "2 Gy x 23 + 1.1 Gy x 20 concomitant boost"),
]
RTOG9003_BEITLER = [
 ("RTOG9003", "RTOG9003", "24613816", "Standard", 2.0, 35, 50, 268, 45.7, False, "median delivered", "Beitler 2014 rate; N from Fu 2000 Table 1"),
 ("RTOG9003", "RTOG9003", "24613816", "Hyperfractionation", 1.2, 68, 50, 263, 53.3, False, "median delivered", "Beitler 2014 rate; N from Fu 2000 Table 1"),
 ("RTOG9003", "RTOG9003", "24613816", "with split", 1.6, 42, 43, 274, 47.8, False, "median delivered", "Beitler 2014 rate; N from Fu 2000 Table 1"),
 ("RTOG9003", "RTOG9003", "24613816", "concomitant boost", 72 / 42, 42, 43, 268, 53.8, False, "median delivered", "Beitler 2014 rate; N from Fu 2000 Table 1"),
]
# Sensitivity sets (decision 2026-10-22: postoperative and NPC only as sensitivity fits).
POSTOP = [
 ("Ang_postop", "Ang_postop", "11597795", "7 weeks", 1.8, 35, 49, 75, 61, False, "planned", "postoperative; 5-y value read from Fig. 3A"),
 ("Ang_postop", "Ang_postop", "11597795", "5 weeks", 1.8, 35, 35, 76, 73, False, "planned", "postoperative; 5-y value read from Fig. 3A"),
 ("Awwad_postop", "Awwad_postop", "11870530", "CF", 2.0, 30, 42, 39, 57, False, "planned", "postoperative"),
 ("Awwad_postop", "Awwad_postop", "11870530", "AHF", 1.4, 33, 12, 31, 88, False, "planned", "postoperative"),
 ("Sanguineti_postop", "Sanguineti_postop", "15708255", "CF", 2.0, 30, 42, 113, 80, False, "median delivered", "postoperative"),
 ("Sanguineti_postop", "Sanguineti_postop", "15708255", "AF", 64 / 35, 35, 35, 113, 78, False, "median delivered", "postoperative; 64 Gy/35 fx"),
 ("pCAIR", "pCAIR", "18342964", "p-CF", 1.8, 35, 48.2, 137, 64, False, "median delivered", "postoperative"),
 ("pCAIR", "pCAIR", "18342964", "p-CAIR", 1.8, 35, 35.6, 137, 70, False, "median delivered", "postoperative"),
]
NPC = [
 ("NPC9902", "NPC9902", "16904519", "CF arm", 2.0, 33, 48, 42, 85, False, "median delivered", "NPC; RT-only pair of the 2x2 factorial; fractions = 66 Gy / 2 Gy"),
 ("NPC9902", "NPC9902", "16904519", "AF arm", 2.0, 33, 41, 52, 78, False, "median delivered", "NPC; RT-only pair of the 2x2 factorial; fractions = 66 Gy / 2 Gy"),
]
# Added 2026-10-24; sensitivity only since 2026-10-25 (crude or first-event competing-risk endpoints).
ADD_2026_10_24 = [
 ("Vienna", "Vienna_CHART", "11054514", "CF", 2.0, 35, 49, 81, 31, False, "planned (7 weeks)", "ultimate crude locoregional control, median follow-up 48 mo"),
 ("Vienna", "Vienna_CHART", "11054514", "V-CHART", 55.3 / 33, 33, 17, 78, 32, False, "planned (17 consecutive days)", "2.5 Gy day 1 then 1.65 Gy twice daily x 32; MMC arm excluded; crude control"),
 ("ELAN_RT", "ELAN_RT", "41713468", "HSC-RT", 2.75, 20, 43, 101, 27.3, True, "median delivered", "frail elderly; split course (2-week planned gap); competing-risk first-event locoregional progression"),
 ("ELAN_RT", "ELAN_RT", "41713468", "SF-RT", 2.0, 35, 52, 100, 24.7, True, "median delivered", "frail elderly; competing-risk first-event locoregional progression"),
 ("Gupta2015", "Gupta2015", "26942141", "Arm-A", 2.0, 33, 38.5, 55, 52.7, False, "planned (5.5 weeks)", "6 fx/week; crude complete response maintained at median follow-up 43 mo"),
 ("Gupta2015", "Gupta2015", "26942141", "Arm-B", 2.0, 33, 45.5, 54, 44.4, False, "planned (6.5 weeks)", "5 fx/week; crude complete response maintained at median follow-up 43 mo"),
]
# Chemoradiotherapy with a common chemotherapy policy in the compared arms (not equal delivered
# exposure: HYPNO planned 4 vs 5 cisplatin cycles): sensitivity fit only.
CHEMO_SAME = [
 ("HYPNO", "HYPNO", "41711624", "HFX", 2.75, 20, 28, 395, 50.7, False, "planned", "cisplatin in 75% of each arm; planned 4 vs 5 weekly cycles"),
 ("HYPNO", "HYPNO", "41711624", "NFX", 2.0, 33, 38.5, 397, 51.2, False, "planned", "cisplatin in 75% of each arm; planned 4 vs 5 weekly cycles"),
 ("ART_DECO", "ART_DECO", "34256319", "ST-IMRT", 65 / 30, 30, 39, 138, 73.4, False, "delivered (as reported)", "same chemotherapy policy, balanced by minimisation"),
 ("ART_DECO", "ART_DECO", "34256319", "DE-IMRT", 67.2 / 28, 28, 37, 138, 71.4, False, "delivered (as reported)", "same chemotherapy policy, balanced by minimisation"),
 ("Shivananjappa", "Shivananjappa", "37006066", "Hypo-SIB", 2.2, 30, 39.4, 25, 92, False, "delivered (as reported)", "cisplatin 100 mg/m2 3-weekly in both arms"),
 ("Shivananjappa", "Shivananjappa", "37006066", "Conv-VMAT", 2.0, 35.5, 50.2, 25, 84, False, "delivered (as reported)", "cisplatin 100 mg/m2 3-weekly in both arms; 35-36 fractions"),
]

# Fraction-size components (number, Gy) for arms that mix fraction sizes within a patient's course.
# Checked in the full texts on 2026-10-25.
COMPONENTS = {
 ("ARTSCAN", "AF (accelerated"): [(23, 2.0), (20, 1.1)],       # 2 Gy + 1.1 Gy per day to 68 Gy
 ("Ghoshal", "Concomitant boost"): [(25, 1.8), (15, 1.5)],     # 45 Gy/25 + 22.5 Gy/15
 ("RTOG8809", "AFX-C"): [(30, 1.8), (11, 1.5)],                # 54 Gy/30 + 1.5 Gy boost on last 11 days
 ("RTOG9003", "concomitant boost"): [(30, 1.8), (12, 1.5)],    # 1.8 Gy/day + 1.5 Gy boost on last 12 days
 ("Vienna", "V-CHART"): [(1, 2.5), (32, 1.65)],                # 2.5 Gy day 1, then 1.65 Gy twice daily
}
# Patient mixtures within an arm: (patients, fractions) per subgroup, same dose per fraction.
MIXTURES = {
 ("KROG0201", "CONV"): [(74, 33), (8, 35)],
 ("KROG0201", "HYPO"): [(65, 28), (9, 30)],
}
# Planned schedule on the common calendar (nominal-clock sensitivity), from the extracted plans.
# For concomitant-boost arms the count is treatment days (boost fractions share a day).
NOMINAL = {
 ("CAIR_2.0Gy", "CAIR, 2.0 Gy"): one(35, EVERY_DAY), ("CAIR_2.0Gy", "Control, 2.0 Gy"): one(35),
 ("CAIR_1.8Gy", "CAIR, 1.8 Gy"): one(40, EVERY_DAY), ("CAIR_1.8Gy", "Control, 1.8 Gy"): one(40),
 ("CAIR_vs_CB", "CAIR"): one(69.1 / 1.8, EVERY_DAY), ("CAIR_vs_CB", "CB"): one(70.1 / 1.8, (1, 2, 1, 1, 2, 0, 0)),
 ("CHART", "CHART"): one(36, (3,) * 7), ("CHART", "Conventional"): one(33),
 ("DAHANCA67", "6 fractions"): one(33.5, MON_SAT), ("DAHANCA67", "5 fractions"): one(33.5),
 ("EORTC22791", "CF"): one(35), ("EORTC22791", "HF"): one(70, BID_MON_FRI),
 ("EORTC22811", "SDF"): one(35),
 # 3 x 1.6 Gy/day for 10 days, 3-4 week rest, boost 3/day for 4-5 days (42 or 45 fractions)
 ("EORTC22811", "MFD (three"): [(1, F, (3, 3, 3, 3, 3, 0, 0), (30, r)) for F in (42, 45) for r in (21, 28)],
 ("EORTC22851", "CF"): one(35),
 # 3 x 1.6 Gy/day: 18 fractions, 12-14 day split, 27 fractions
 ("EORTC22851", "AF"): [(1, 45, (3, 3, 3, 3, 3, 0, 0), (18, r)) for r in (12, 14)],
 ("Ghoshal", "Concomitant boost"): one(25), ("Ghoshal", "Conventional"): one(33),   # 25 treatment days
 ("Hliniak", "CF"): one(33), ("Hliniak", "AF"): one(33, (1, 1, 1, 2, 1, 0, 0)),
 ("IAEA_ACC", "six fractions"): one(34, MON_SAT), ("IAEA_ACC", "five fractions"): one(34),
 ("Jackson", "Once daily"): one(33), ("Jackson", "Twice daily"): one(33, BID_MON_FRI),
 ("KROG0201", "CONV"): [(74, 33, MON_FRI, None), (8, 35, MON_FRI, None)],
 ("KROG0201", "HYPO"): [(65, 28, MON_FRI, None), (9, 30, MON_FRI, None)],
 ("PMH_HF", "HF"): one(40, BID_MON_FRI), ("PMH_HF", "CF"): one(20),
 ("RTOG8313_A", "67.2"): one(56, BID_MON_FRI), ("RTOG8313_A", "72.0"): one(60, BID_MON_FRI),
 ("RTOG8313_A", "76.8"): one(64, BID_MON_FRI),
 ("RTOG8809", "AHFX-S"): one(42, BID_MON_FRI, (24, 14)), ("RTOG8809", "AFX-C"): one(30),
 ("RTOG9003", "Standard"): one(35), ("RTOG9003", "Hyperfractionation"): one(68, BID_MON_FRI),
 ("RTOG9003", "with split"): one(42, BID_MON_FRI, (24, 14)), ("RTOG9003", "concomitant boost"): one(30),
 ("RTOG9512", "SFX"): one(35), ("RTOG9512", "HFX"): one(66, BID_MON_FRI),
 ("TROG9101", "CRT"): one(35), ("TROG9101", "ART"): one(33, (2, 2, 2, 2, 1, 0, 0)),
 ("Yamazaki_small", "Arm A-1"): one(30), ("Yamazaki_small", "Arm B-1"): one(25),
 ("Yamazaki_large", "Arm A-2"): one(33), ("Yamazaki_large", "Arm B-2"): one(28),
 ("ARTSCAN", "CF (conventional"): one(34), ("ARTSCAN", "AF (accelerated"): one(23),
 ("Vienna", "CF"): one(35), ("Vienna", "V-CHART"): one(17, EVERY_DAY),   # 17 consecutive treatment days (33 fractions)
 ("ELAN_RT", "HSC-RT"): one(20, MON_FRI, (10, 14)), ("ELAN_RT", "SF-RT"): one(35),
 ("Gupta2015", "Arm-A"): one(33, MON_SAT), ("Gupta2015", "Arm-B"): one(33),
 ("HYPNO", "HFX"): one(20), ("HYPNO", "NFX"): one(33, MON_SAT),
 ("ART_DECO", "ST-IMRT"): one(30), ("ART_DECO", "DE-IMRT"): one(28),
 ("Shivananjappa", "Hypo-SIB"): one(30), ("Shivananjappa", "Conv-VMAT"): one(35.5),
}
# Rows whose N is not the matched extraction row's n_analyzed, with the reason.
N_SOURCE = {
 ("RTOG9003", "24613816"): "Beitler 2014 gives no per-arm N (1076 analyzed); Fu 2000 Table 1 N used (sum 1073)",
}


def build(argv, matches=None):
    """matches: optional list; receives (output row, matched extraction row) pairs for provenance."""
    R = json.load(open(os.path.join(EXTRACT, "all_extract.json"), encoding="utf-8"))
    R = R + json.load(open(os.path.join(EXTRACT, "supp_extract_2026_10_24.json"), encoding="utf-8"))
    arms = list(A)
    if "--rtog9003-beitler" in argv:
        arms = [a for a in arms if a[0] != "RTOG9003"] + RTOG9003_BEITLER
    if "--artscan" in argv:
        arms = arms + ARTSCAN
    if "--add-2026-10-24" in argv:
        arms = arms + ADD_2026_10_24
    if "--chemo" in argv:
        arms = arms + CHEMO_SAME
    if "--drop-elan" in argv:
        arms = [a for a in arms if a[0] != "ELAN_RT"]
    if "--no-figure" in argv:
        arms = [a for a in arms if a[0] not in ("CHART", "EORTC22811")]
    if "--postop" in argv:
        arms = arms + POSTOP
    if "--npc" in argv:
        arms = arms + NPC
    nscale = 1.0
    if "--n-scale" in argv:
        nscale = float(argv[argv.index("--n-scale") + 1])

    out, problems = [], []
    for unit, trial, pmid, arm, d, nfx, T, N, v, fail, trule, note in arms:
        m = [r for r in R if r["trial"] == trial and str(r["pmid"]) == pmid and arm in str(r["arm"])
             and r.get("value_pct") not in ("", None) and abs(float(r["value_pct"]) - v) < 1e-9]
        if not m:
            problems.append(f"{unit} {arm}: no extracted row of {pmid} with value {v}")
            continue
        # Several matching extraction rows must agree on everything used or kept (review R08).
        sig = lambda x: tuple(str(x.get(k, "")) for k in ("arm", "n_analyzed", "value_ci", "endpoint", "estimate_type",
                                                          "timepoint_years", "value_source"))
        if len({sig(x) for x in m}) > 1:
            problems.append(f"{unit} {arm}: {len(m)} extraction rows of {pmid} match with different content")
        r = m[0]
        qok = r.get("_quote_ok") or {}
        if r.get("value_source") not in ("figure_read", "table_pdf_only"):
            if not qok:
                problems.append(f"{unit} {arm}: source row has no recorded quote check")
            elif not all(qok.values()):
                problems.append(f"{unit} {arm}: source row failed the quote check")
        # Patient number must match the report the rate comes from, or carry a stated reason.
        nsrc = "n_analyzed of the matched extraction row"
        try:
            n_ex = int(float(r.get("n_analyzed")))
        except (TypeError, ValueError):
            n_ex = None
        if n_ex != N:
            if (trial, pmid) in N_SOURCE:
                nsrc = N_SOURCE[(trial, pmid)]
            else:
                problems.append(f"{unit} {arm}: N {N} differs from extraction n_analyzed {r.get('n_analyzed')!r}")
        key = next((k for k in COMPONENTS if k[0] == unit and k[1] in arm), None)
        comps = COMPONENTS[key] if key else [(nfx, d)]
        if abs(sum(c[0] for c in comps) - nfx) > 1e-9 or abs(sum(c[0] * c[1] for c in comps) - nfx * d) > 1e-6:
            problems.append(f"{unit} {arm}: fraction components {comps} do not sum to {nfx} x {d}")
        mkey = next((k for k in MIXTURES if k[0] == unit and k[1] == arm), None)
        if mkey:
            mix = MIXTURES[mkey]
            if sum(p for p, _ in mix) != N or abs(sum(p * f for p, f in mix) / N - nfx) > 1e-9:
                problems.append(f"{unit} {arm}: patient mixture {mix} does not give N {N} and {nfx} fractions")
        nkey = next((k for k in NOMINAL if k[0] == unit and k[1] == arm), None)
        if nkey is None and unit not in {a[0] for a in POSTOP + NPC}:
            problems.append(f"{unit} {arm}: no planned calendar for the nominal-clock set")
        Tnom = nominal(NOMINAL[nkey]) if nkey else ""
        tcp = 100 - v if fail else v
        n_out = N if nscale == 1.0 else max(1, int(round(N * nscale)))
        use_nom = "--nominal-clock" in argv and Tnom != ""
        Tuse = Tnom if use_nom else T
        out.append(dict(
            row=len(out) + 1, first_author=trial, ref_no=unit, arm=r["arm"], dose_per_fraction_gy=round(d, 4),
            n_fractions=round(nfx, 3), total_time_days=Tuse, n_patients=n_out, tcp_actual_pct=round(tcp, 2),
            endpoint=r["endpoint"] + (" (converted: 100 - failure)" if fail else ""),
            timepoint_years=r["timepoint_years"], value_source=r["value_source"], pmid=pmid,
            location=r.get("location", ""),
            time_rule="nominal calendar (Monday start, inclusive)" if use_nom else trule, notes=note,
            total_dose_gy=round(sum(c[0] * c[1] for c in comps), 6),
            sum_d2_gy2=round(sum(c[0] * c[1] ** 2 for c in comps), 6),
            fraction_components=";".join(f"{c[0]:g}x{c[1]:g}" for c in comps) if key else "",
            T_source_days=T, T_source_rule=trule, T_nominal_days=Tnom,
            estimate_type=r.get("estimate_type", ""), endpoint_definition=r.get("endpoint_definition", ""),
            source_value_pct=v, source_value_is_failure=fail, n_source=nsrc,
            n_scale=nscale, n_reported=N))
        if matches is not None:
            matches.append((out[-1], r))
    return out, problems


def km_only(out):
    """Keep only fitting units whose arms are all Kaplan-Meier estimates; renumber rows."""
    bad = {o["ref_no"] for o in out if o["estimate_type"] != "actuarial/KM"}
    out = [o for o in out if o["ref_no"] not in bad]
    for i, o in enumerate(out, 1):
        o["row"] = i
    return out


def main(argv):
    out, problems = build(argv)
    if problems:
        sys.exit("REFUSED:\n" + "\n".join(problems))
    if "--km-only" in argv:
        out = km_only(out)
    with open(argv[1], "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
        w.writeheader(); w.writerows(out)
    print(len(out), "arms,", len({o["ref_no"] for o in out}), "units,",
          len({o["first_author"] for o in out}), "trials ->", argv[1])


if __name__ == "__main__":
    main(sys.argv)
