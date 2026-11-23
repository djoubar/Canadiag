# ==============================================================================
#
#        MIMIC-IV : 1 LIGNE PAR HÉMOCULTURE, CANDIDA POSITIVE / NEGATIVE
#        + covariables mesurées AVANT l'hémoculture (même hospitalisation)
#
# ==============================================================================
#
# Unité d'analyse : un prélèvement d'hémoculture (micro_specimen_id)
# Outcome         : "Positive" si au moins un isolat Candida / levure, sinon "Negative"
#
# Fenêtre des covariables : [début d'hospitalisation (edregtime ou admittime) ;
#                            heure de l'hémoculture[
#   - chirurgie         : passage en PACU (horodaté) OU acte CCAM-like "bloc"
#                         (ICD-10-PCS section 0, abord ouvert/cœlioscopique ;
#                          ICD-9 01-86 hors gestes de chevet/endoscopies)
#   - chirurgie_abdo    : idem restreint au digestif / péritoine
#   - transfusion       : inputevents (CGR, PFC, plaquettes, cryo, sang total,
#                         OR/PACU) OU code ICD transfusion
#   - creat_max         : créatinine sanguine max (mg/dL), labevents 50912/52546
#   - eer               : procedureevents (CRRT, CVVH(D)(F), HD, SCUF, DP),
#                         inputevents spécifiques CRRT, OU code ICD dialyse
#
# Les actes ICD n'ont qu'une date (pas d'heure) : ils ne comptent que s'ils
# sont datés d'un jour STRICTEMENT antérieur à l'hémoculture (pas de fuite).
# ==============================================================================

library(data.table)

mimic_dir <- "mimic_iv" # racine contenant core/, hosp/, icu/
out_dir <- file.path(mimic_dir, "derived")
dir.create(out_dir, showWarnings = FALSE)

rd <- function(path, ...) fread(cmd = paste("zcat", file.path(mimic_dir, path)), ...)

# ------------------------------------------------------------------------------
# 1. Hémocultures
# ------------------------------------------------------------------------------

micro <- rd(
  "hosp/microbiologyevents.csv.gz",
  select = c(
    "subject_id",
    "hadm_id",
    "micro_specimen_id",
    "chartdate",
    "charttime",
    "spec_type_desc",
    "org_name"
  )
)

bc_types <- c("BLOOD CULTURE", "BLOOD CULTURE ( MYCO/F LYTIC BOTTLE)")
candida_regex <- "CANDIDA|^YEAST" # en hémoculture, "YEAST" = Candida non encore identifié

bc_raw <- micro[spec_type_desc %in% bc_types]

bc <- bc_raw[,
  .(
    subject_id = subject_id[1],
    hadm_id = hadm_id[1],
    spec_type_desc = spec_type_desc[1],
    # sans heure -> minuit du jour (conservateur : aucune covariable du jour même)
    t_bc = if (!is.na(charttime[1])) charttime[1] else chartdate[1],
    candida_org = paste(unique(org_name[grepl(candida_regex, org_name)]), collapse = "; "),
    cancelled = all(org_name == "CANCELLED")
  ),
  by = micro_specimen_id
]

n_cancel <- bc[cancelled == TRUE, .N]
bc <- bc[cancelled == FALSE][, cancelled := NULL]
bc[,
  candida := factor(
    ifelse(candida_org != "", "Positive", "Negative"),
    levels = c("Negative", "Positive")
  )
]
bc[candida_org == "", candida_org := NA_character_]

# ------------------------------------------------------------------------------
# 2. Rattachement à l'hospitalisation (hadm_id manquant = souvent prélevé aux urgences)
# ------------------------------------------------------------------------------

adm <- rd(
  "core/admissions.csv.gz",
  select = c("subject_id", "hadm_id", "admittime", "dischtime", "edregtime")
)
adm[, adm_start := pmin(admittime, edregtime, na.rm = TRUE)]

miss <- bc[is.na(hadm_id), .(micro_specimen_id, subject_id, t_bc)]
found <- adm[
  miss,
  on = .(subject_id, adm_start <= t_bc, dischtime >= t_bc),
  .(micro_specimen_id = i.micro_specimen_id, hadm_id = x.hadm_id),
  nomatch = NULL,
  mult = "first"
]
bc[found, on = "micro_specimen_id", hadm_id := i.hadm_id]

bc[adm, on = "hadm_id", `:=`(win_start = i.adm_start, dischtime = i.dischtime)]

# Séjour en réanimation au moment du prélèvement
icu <- rd("icu/icustays.csv.gz", select = c("subject_id", "stay_id", "intime", "outtime"))
icu_hit <- icu[
  bc,
  on = .(subject_id, intime <= t_bc, outtime >= t_bc),
  .(micro_specimen_id = i.micro_specimen_id, stay_id = x.stay_id, icu_intime = x.intime),
  nomatch = NULL,
  mult = "first"
]
bc[icu_hit, on = "micro_specimen_id", `:=`(stay_id = i.stay_id, icu_intime = i.icu_intime)]
bc[, icu_h_before_bc := round(as.numeric(difftime(t_bc, icu_intime, units = "hours")), 1)]

# ------------------------------------------------------------------------------
# 3. Helpers : évènement dans [win_start ; t_bc[
# ------------------------------------------------------------------------------

keys <- bc[!is.na(win_start), .(micro_specimen_id, subject_id, win_start, t_bc)]

# ev : data.table(subject_id, t) -> vecteur logique aligné sur bc
any_before <- function(ev) {
  hit <- unique(ev[
    keys,
    on = .(subject_id, t >= win_start, t < t_bc),
    .(micro_specimen_id = i.micro_specimen_id),
    nomatch = NULL
  ])
  out <- bc$micro_specimen_id %in% hit$micro_specimen_id
  out[is.na(bc$win_start)] <- NA # hémoculture hors hospitalisation
  as.integer(out)
}

# Actes ICD : date seule -> fin de journée, donc compte seulement si jour < jour du prélèvement
icd_ev <- function(codes_dt) {
  codes_dt[, .(subject_id, t = as.POSIXct(chartdate, tz = "UTC") + 86399)]
}

# ------------------------------------------------------------------------------
# 4. Sources
# ------------------------------------------------------------------------------

pr_icd <- rd(
  "hosp/procedures_icd.csv.gz",
  select = c("subject_id", "chartdate", "icd_code", "icd_version"),
  colClasses = list(character = "icd_code")
)
pr9 <- pr_icd[icd_version == 9]
pr10 <- pr_icd[icd_version == 10]

# --- Chirurgie ----------------------------------------------------------------
# ICD-9 : chapitres opératoires 01-86, hors gestes de chevet / diagnostiques
icd9_non_or <- c(
  "3891",
  "3892",
  "3893",
  "3894",
  "3895",
  "3897",
  "3898",
  "3899", # cathéters / ponctions vasculaires
  "3995",
  "5498", # dialyse
  "3321",
  "3322",
  "3323",
  "3324",
  "3327", # bronchoscopies
  "4213",
  "4219",
  "4222",
  "4223",
  "4224",
  "4411",
  "4413",
  "4414", # endoscopies hautes
  "4513",
  "4514",
  "4516",
  "4521",
  "4522",
  "4523",
  "4524",
  "4525", # endoscopies basses
  "5110",
  "5111",
  "5114",
  "5210",
  "5213",
  "5114", # CPRE
  "5491",
  "3491",
  "3404",
  "0331",
  "8622",
  "8628",
  "8604",
  "8605" # ponctions, drains, parage
)
is_or9 <- function(code) {
  ch <- suppressWarnings(as.integer(substr(code, 1, 2)))
  !is.na(ch) & ch >= 1 & ch <= 86 & !(substr(code, 1, 4) %in% icd9_non_or)
}
# ICD-10-PCS : section 0 (Medical & Surgical), abord ouvert (0) ou endoscopique percutané (4)
is_or10 <- function(code) substr(code, 1, 1) == "0" & substr(code, 5, 5) %in% c("0", "4")

surg_icd <- rbind(pr9[is_or9(icd_code)], pr10[is_or10(icd_code)])

abdo9 <- function(code) {
  ch <- suppressWarnings(as.integer(substr(code, 1, 2)))
  !is.na(ch) & ch >= 42 & ch <= 54
}
abdo10 <- function(code) {
  substr(code, 2, 2) %in%
    c("D", "F") | # digestif, hépato-bilio-pancréatique
    (substr(code, 2, 2) == "W" & substr(code, 4, 4) %in% c("F", "G", "H", "J", "P")) # paroi, péritoine, rétropéritoine
}
abdo_icd <- surg_icd[(icd_version == 9 & abdo9(icd_code)) | (icd_version == 10 & abdo10(icd_code))]

transfers <- rd("core/transfers.csv.gz", select = c("subject_id", "careunit", "intime"))
pacu <- transfers[careunit == "PACU", .(subject_id, t = intime)]

bc[, surgery := any_before(rbind(pacu, icd_ev(surg_icd)))]
bc[, surgery_abdo := any_before(icd_ev(abdo_icd))]

# --- Transfusion --------------------------------------------------------------
blood_items <- c(
  225168,
  220996,
  226368,
  227070, # CGR (réa, OR, PACU)
  220970,
  220971, # PFC
  225170,
  226369,
  227071, # plaquettes
  225171,
  226371, # cryoprécipités
  221013 # sang total
)
crrt_input_items <- c(227525, 227536, 230044) # Ca, KCl, héparine "CRRT"

inputs <- rd("icu/inputevents.csv.gz", select = c("subject_id", "starttime", "itemid"))
inputs <- inputs[itemid %in% c(blood_items, crrt_input_items)]

transf_icd <- rbind(
  pr9[startsWith(icd_code, "990")], # 99.0x transfusions
  pr10[startsWith(icd_code, "302")] # 302xxxx transfusion (circulatoire)
)
bc[,
  transfusion := any_before(rbind(
    inputs[itemid %in% blood_items, .(subject_id, t = starttime)],
    icd_ev(transf_icd)
  ))
]

# --- Épuration extra-rénale ---------------------------------------------------
rrt_proc_items <- c(225802, 225803, 225809, 225955, 225441, 225805)
proc <- rd("icu/procedureevents.csv.gz", select = c("subject_id", "starttime", "itemid"))

eer_icd <- rbind(
  pr9[icd_code %in% c("3995", "5498")],
  pr10[startsWith(icd_code, "5A1D") | icd_code == "3E1M39Z"]
)
bc[,
  eer := any_before(rbind(
    proc[itemid %in% rrt_proc_items, .(subject_id, t = starttime)],
    inputs[itemid %in% crrt_input_items, .(subject_id, t = starttime)],
    icd_ev(eer_icd)
  ))
]

# --- Créatinine max -----------------------------------------------------------
# préfiltre awk : labevents complet = plusieurs Go
crea <- fread(
  cmd = sprintf(
    "zcat %s | awk -F, 'NR==1 || $5==50912 || $5==52546'",
    file.path(mimic_dir, "hosp/labevents.csv.gz")
  ),
  select = c("subject_id", "charttime", "valuenum")
)[!is.na(valuenum)]

crea_max <- crea[
  keys,
  on = .(subject_id, charttime >= win_start, charttime < t_bc),
  .(
    micro_specimen_id = i.micro_specimen_id,
    creat_max = if (all(is.na(x.valuenum))) NA_real_ else max(x.valuenum, na.rm = TRUE)
  ),
  by = .EACHI
][, .(micro_specimen_id, creat_max)]
bc[crea_max, on = "micro_specimen_id", creat_max := i.creat_max]

# ------------------------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------------------------

df_bc <- bc[, .(
  micro_specimen_id,
  subject_id,
  hadm_id,
  stay_id,
  spec_type_desc,
  t_bc,
  icu_h_before_bc,
  candida,
  candida_org,
  surgery,
  surgery_abdo,
  transfusion,
  creat_max,
  eer
)][order(subject_id, t_bc)]

saveRDS(as.data.frame(df_bc), file.path(out_dir, "bc_candida.rds"))
fwrite(df_bc, file.path(out_dir, "bc_candida.csv"))

message(sprintf(
  "%d hémocultures (%d annulées exclues), %d patients, %d Candida positives ; %d hors hospitalisation (covariables NA)",
  nrow(df_bc),
  n_cancel,
  uniqueN(df_bc$subject_id),
  sum(df_bc$candida == "Positive"),
  sum(is.na(df_bc$hadm_id))
))
