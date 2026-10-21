##LPO IV extraction and MR
##Install required packages
install.packages("data.table")
install.packages("TwoSampleMR")
install.packages("remotes")
remotes::install_github("MRCIEU/genetics.binaRies")
library(TwoSampleMR)
library(data.table)
library(genetics.binaRies)
library(ieugwasr)

##loading dataset
lpoinstruments <- fread(
  "/Users/laurabell/Documents/workdir/lpoinstruments.tsv"
)
names(lpoinstruments)
head(lpoinstruments)

##filtering significant instruments
LPOsig <- lpoinstruments[p_value < 5e-8]
nrow(LPOsig)
unique(LPOsig$chromosome)

##defining LPO gene region
lpo_start <- 58218548
lpo_end <- 58268518

##defining cis and trans SNPs
lposigcis <- LPOsig[
  chromosome == 17 &
    base_pair_location >= (lpo_start - 500000) &
    base_pair_location <= (lpo_end + 500000)
]

nrow(lposigcis)
head(lposigcis)

##rename columns
lposigcis_exposure <- lposigcis %>%
  filter(!is.na(rs_id)) %>%
  transmute(
    SNP = rs_id,
    beta = beta,
    se = standard_error,
    effect_allele = effect_allele,
    other_allele = other_allele,
    eaf = effect_allele_frequency,
    pval = p_value,
    chr = chromosome,
    pos = base_pair_location,
    samplesize = n
  )

##creating data frame
lposigcis_exposure <- as.data.frame(lposigcis_exposure)

lposigcis_exposure <- format_data(
  lposigcis_exposure,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "pos",
  samplesize_col = "samplesize"
)

head(lposigcis_exposure)
lposigcis_exposure$exposure <- "LPO"
lposigcis_exposure$id.exposure <- "LPO"

lposigcis_clumped <- clump_data(
  lposigcis_exposure,
  clump_kb = 500,
  clump_r2 = 0.01,
  clump_p1 = 5e-8,
  pop = "EUR"
)

##clumping with plink

plink <- genetics.binaRies::get_plink_binary()
plink
http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
ld_reference <- "/Users/laurabell/Documents/workdir/1kg.v3/EUR"

head(lposigcis)
nrow(lposigcis_exposure)

clump_input <- data.frame(
  rsid = lposigcis_exposure$SNP,
  pval = lposigcis_exposure$pval.exposure,
  id = "LPO"
)
lpo_clumped <- ieugwasr::ld_clump(
  dat = clump_input,
  clump_kb = 500,
  clump_r2 = 0.01,
  clump_p = 5e-8,
  bfile = ld_reference,
  plink_bin = plink
)

head(lpo_clumped)
names(lpo_clumped)

lposigcis_clumped

lpo_rsid <- lpo_clumped$rsid
lpo_rsid

##selecting clumped instruments from exposure dataset

lpo_instruments_filtered <- lposigcis_exposure %>%
  dplyr::filter(
    SNP == lpo_rsid,
  )

head(lpo_instruments_filtered)


##extracting outcomes

crc_outcome <- extract_outcome_data(
  snps = lpo_instruments_filtered$SNP,
  outcomes = "ieu-b-4965"
)
head(lpo_crc)

##harmonising dataset

lpo_crc <- harmonise_data(
  exposure_dat = lpo_instruments_filtered,
  outcome_dat = crc_outcome,
  action = 2)

head(lpo_crc)  

##doing MR

lpo_mr <- mr(lpo_crc_harmonised)
lpo_mr

##making odds ratio
lpo_odds_ratios<-generate_odds_ratios(lpo_mr)
write.csv(lpo_odds_ratios, "lpo_odds_ratios.csv", row.names = FALSE)

##creating graphs for LPO
lpo_dat <- subset(
  lpo_mr,
  exposure == "LPO"
)
lpo_dat <- lpo_dat[!duplicated(lpo_dat$SNP), ]
lpo_dat[, c(
  "SNP",
  "beta.exposure",
  "beta.outcome",
  "se.exposure",
  "se.outcome"
)]

mr_scatter_plot(lpo_mr, lpo_crc)
res_singlelpo <- mr_singlesnp(lpo_crc)
mr_forest_plot(res_singlelpo)
mr_funnel_plot(res_singlelpo)

unique(lpo_mr$method)
dim(lpo_dat)
head(lpo_dat)

nrow(lpo_dat)

lpo_mr

##calculating wald ratios
lpo_crc$wald_ratio <- lpo_crc$beta.outcome / lpo_crc$beta.exposure
lpo_crc$wald_ratio


