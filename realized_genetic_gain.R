library(tidyverse)
library(asreml)

################################################################################
################################################################################
#### read the data and qa/qc (remove empty reps and constant values)
################################################################################
################################################################################

#### read the data
dat <- read_csv("data/cassava_uyt_2016_2021.csv")
str(dat)
unique(dat$yearTesting)
length(unique(dat$studyName))
length(unique(dat$loc))
length(unique(dat$geno))
unique(dat$design)


#### remove empty reps and constant values
datCleaned <- dat |>
  mutate(across(fyld, ~ replace(.x, is.nan(.x), NA)))  |>
  nest_by(studyName, rep) |>
  mutate(
    missing = list(sum(is.na(data$fyld))),
    missing = round(100 * missing / nrow(data), 0),
    variance = var(data$fyld, na.rm = TRUE)
  ) |>
  filter(missing != 100) |>
  filter(!is.na(variance)) |>
  filter(variance != 0) |>
  filter(missing < 30) |>
  select(-c(missing, variance)) |>
  unnest(data) |>
  ungroup() |>
  nest_by(studyName) |>
  mutate(nRep = length(unique(data$rep))) |>
  filter(nRep > 1) |>
  select(-nRep) |>
  unnest(data) |>
  ungroup() |>
  mutate(across(
    c(studyName, rep, block, yearTesting, trial, loc, geno, entryType),
    as.factor
  ))

##################################################################################################
##################################################################################################
#### STEP 1: estimate and discard env. with low heritability (Cullis)
#### Linear Mixed Model (LMM) for each env.
##################################################################################################
##################################################################################################

#### Run single trial analysis: geno as random
single.study.random <- datCleaned |>
  nest_by(studyName, design, yearTesting) |>
  mutate(
    data = list(droplevels(data)),
    model.rando = ifelse(design == "Alpha",
                         list(tryCatch(
                           asreml(
                             fyld ~ rep,
                             random = ~ geno + block,
                             workspace = "3gb",
                             data = data
                           ),
                           error = function(e) {
                             NULL
                           }
                         )),
                         list(tryCatch(
                           # RCBD
                           asreml(
                             fyld ~ rep,
                             random = ~ geno,
                             workspace = "3gb",
                             data = data
                           ),
                           error = function(e) {
                             NULL
                           }
                         )))
  )

#### Filter null and estimate H2
single.study.random <- single.study.random |>
  filter(!is.null(model.rando)) |>
  mutate(
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    vg = list(
      summary(model.rando)$varcomp |>
        rownames_to_column(var = "Effect") |>
        select(Effect, everything()) |>
        filter(grepl("geno", Effect)) |>
        pull(component)
    ),
    vd.mat = list(
      predict(
        model.rando,
        classify = "geno",
        only = "geno",
        sed = TRUE
      )$sed ^ 2
    ),
    vd.avg = list(mean(vd.mat[upper.tri(vd.mat, diag = FALSE)])),
    H2 = list(1 - (vd.avg / (vg * 2)))
  )

##################################################################################################
##################################################################################################
#### STEP 2: example of a one-stage approach
#### Run a combined Linear Mixed Model (LMM) and estimate combined BLUEs and weights
##################################################################################################
##################################################################################################

#### Filter low H2 and prepare the data
datCombined <- single.study.random |>
  filter(H2 >= 0.2) |>
  select(studyName, design, yearTesting, data) |>
  unnest(data) |>
  ungroup()
datCombined <- droplevels(datCombined)

#### Run combined trial analysis: geno as fixed
model.combined <- tryCatch(
  asreml(
    fixed = fyld ~ yearTesting + geno,
    random = ~ loc + studyName + geno:loc + geno:yearTesting + geno:studyName,
    residual = ~ dsum( ~ units | studyName),
    workspace = "3gb",
    data = datCombined
  ),
  error = function(e) {
    NULL
  }
)

if (!is.null(model.combined)) {
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  
  #### estimate combined BLUEs and weights
  predictions.combined <- predict.asreml(model.combined, pworkspace = "5gb", classify = "geno")
  blues.combined <- predictions.combined$pvals |>
    mutate(weight = 1 / (std.error) ^ 2) |>
    left_join(
      dat |>
        select(geno, yearCloning) |>
        filter(yearCloning != 0) |>
        droplevels() |>
        group_by(geno) |>
        summarize(yearCloning = min(yearCloning, na.rm = TRUE))
    ) |>
    filter(predicted.value > 0) |>
    filter(!is.na(yearCloning)) |>
    filter(!is.na(predicted.value)) |>
    left_join(dat |>
                select(geno, entryType) |>
                distinct(geno, .keep_all = TRUE)) |>
    filter(entryType != "check") |>
    droplevels()
  
  
  ##################################################################################################
  ##################################################################################################
  #### STEP 3: weighted-linear regression of the combined BLUEs on the year of origin
  ##################################################################################################
  ##################################################################################################
  
  nY = length(unique(blues.combined$yearCloning))
  
  if (nY >= 5) {
    fit.regression <- lm(predicted.value ~ yearCloning,
                         weights = weight,
                         data = blues.combined)
    slope <-  fit.regression$coefficients[2]
    intercept <-  fit.regression$coefficients[1]
    first.Year.geno <-  min(blues.combined$yearCloning)
    last.Year.geno <-  max(blues.combined$yearCloning)
    baseline <-  intercept + first.Year.geno * slope
    gg.percent <-  100 * slope / baseline
    pvalue <- summary(fit.regression)$coefficients |>
      as.data.frame() |>
      rownames_to_column("Effect") |>
      filter(Effect == "yearCloning") |>
      pull("Pr(>|t|)")
    slope.sd <- summary(fit.regression)$coefficients |>
      as.data.frame() |>
      rownames_to_column("Effect") |>
      filter(Effect == "yearCloning") |>
      pull("Std. Error")
    
    out.gg <- tibble(
      trial = "UYT",
      trait = "fyld",
      slope = slope,
      slope.sd = slope.sd,
      pvalue = pvalue,
      gg.percent = gg.percent
    )
  }
  
  out.gg
  
}
