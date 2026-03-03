
#####################Step 1: Data Cleaning###################

library(chronosphere) 
library(divDyn) 
library(fossilbrush) 
library(tidyverse)

#fetch pbdb data, filter only neccessary things, apply marine filter

pbdb <- fetch("pbdb")
attributes(pbdb)$chronosphere$API

sel <- c("occurrence_no", "collection_no", "collection_name", "cc", "identified_rank", "identified_name",  
         "accepted_rank", "accepted_name", 
         "early_interval", "late_interval", "max_ma", "min_ma", "reference_no", 
         "phylum", "class", "order", "family", "genus", 
         "lng", "lat", "cc","motility", "life_habit", "diet",
         "abund_value", "environment")

dat <- subset(pbdb, select=sel)
setwd("C:\\Users\\fwiss\\Documents\\Uni\\BA_material")
source("FilterMarine.R")

dat <- subset(dat, accepted_rank!="")
dat <- subset(dat, genus!="") 
dat <- subset(dat, life_habit!="amphibious")

data(keys)

stgMin <- categorize(dat$early_interval, keys$stgInt)
stgMax <- categorize(dat$late_interval, keys$stgInt)

stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

dat$stg <- NA
stgCondition <- c(which(stgMax == stgMin), which(stgMax == -1))
dat$stg[stgCondition] <- stgMin[stgCondition]

#fossilbrush, correct typos and synonyms

dat <- dat %>%
  chrono_scale(
    srt = "early_interval", 
    end = "late_interval", 
    max_ma = "max_ma",      
    min_ma = "min_ma",      
    verbose = FALSE         
  ) %>%
  filter(max_ma > min_ma)

corrections <- list(
  "Hipparionix" = "Hipparionyx",
  "Sacamminidae" = "Saccamminidae",
  "Globorotalitidae" = "Globorotaliidae",
  "Paleonisciformes" = "Palaeonisciformes",
  "Ichthyolariidae" = "Ichtyolariidae",     
  "Callorhinchidae" = "Callorhynchidae",
  "Dermochelyoidae" = "Dermochelyidae",
  "Acrocrinitidae" = "Acrocrinidae",
  "Stropheodontidae" = "Strophodontidae",
  "Palmatolepididae" = "Palmatolepidae",
  "Helodontiformes" = "Hybodontiformes",
  "Cassiduloidea" = "Cassiduloida",         
  "Pholadida" = "Pholadomyida",
  "Stromatoporellida" = "Stromatoporida",
  "Craniopsida" = "Craniida",
  "Conularina" = "Conulariida",
  "Productinidae" = "Productidae",
  "Zaphrentoididae" = "Zaphrentidae",
  "Osteolepiformes" = "Osteolepidiformes",
  "Meganterididae" = "Megathyrididae",      
  "Nerinellidae" = "Nerineidae",
  "Opisthonematidae" = "Orthonematidae",
  "Helicocystitidae" = "Hemicystitidae"
)

clean_taxa_names <- function(x) {
  if (all(is.na(x))) return(x)
  x <- as.character(x)
  x <- str_remove(x, "\\s*\\(.*\\).*") 
  x <- str_remove_all(x, "(?i)[?†]") 
  x <- str_remove_all(x, "(?i)\\b(cf\\.|aff\\.|n\\. gen\\.|gen\\.|subgen\\.)\\s*")
  x <- str_remove_all(x, "(?i)\\s+(sp\\.|spp\\.|indet\\.|subsp\\.|var\\.)$")
  x <- str_remove_all(x, "[^a-zA-Z\\-]")
  x <- str_trim(x)
  x[x == ""] <- NA
  return(x)
}

my_ranks <- c("phylum", "class", "order", "family", "genus")

dat <- dat %>%
  mutate(across(all_of(my_ranks), clean_taxa_names))

for (rank in my_ranks) {
  matches <- dat[[rank]] %in% names(corrections)
  if (any(matches)) {
    print(paste("correct to", rank, ":", unique(dat[[rank]][matches])))
    dat[[rank]][matches] <- unlist(corrections[dat[[rank]][matches]])
  }
}

bad_placeholders <- c("NO_GENUS_SPECIFIED", "NO_FAMILY_SPECIFIED", 
                      "NO_ORDER_SPECIFIED", "NO_CLASS_SPECIFIED", 
                      "Indet", "Indeterminate")

dat <- dat %>%
  mutate(across(all_of(my_ranks), ~ replace(., . %in% bad_placeholders, NA)))

my_suffixes <- list(NULL, NULL, NULL, NULL, c("ina", "ella", "etta"))

check_final <- check_taxonomy(
  dat, 
  suff_set = my_suffixes, 
  ranks = my_ranks,
  clean_name = TRUE,       
  resolve_duplicates = TRUE, 
  jump = 5,                
  verbose = TRUE           
)

dat_clean <- check_final$data
print(summary(check_final))

dat <- dat_clean

dens_matrix <- densify(dat)

pacm_res <- pacmacro_ranges(
  dat, 
  tail.flag = 0.35, 
  rank = "genus", 
  srt = "max_ma", 
  end = "min_ma"
)

pranges <- pacm_res$kdensity
suspect_taxa <- pranges[pranges$tflag0.35 == 1, ]

if (nrow(suspect_taxa) > 0) {
  correction_table <- data.frame(
    genus = rownames(suspect_taxa),
    max_ma = suspect_taxa$FAD95,
    min_ma = suspect_taxa$LAD95
  )
  
  pflags <- flag_ranges(dat, correction_table, verbose = FALSE) 
  
  dat$strat_error_code <- pflags$occurrence_status
  
} else {
  dat$strat_error_code <- "OK"
}

itp <- threshold_ranges(
  dat, 
  win = 10,   
  thresh = 10,    
  rank = "genus", 
  srt = "max_ma", 
  end = "min_ma"
)

dat$suggested_split_taxon <- itp$data

#merging with body size data

heim_data <- read.delim("C:\\Users\\fwiss\\Documents\\Uni\\Bachelor_Analyse\\Heim2015_supplementary_data_file.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(heim_data)

if("taxon_name" %in% colnames(heim_data)) {
  heim_prep <- heim_data %>% rename(genus = taxon_name)
} else {
  heim_prep <- heim_data 
}


#median size per genus

heim_sizes <- heim_prep %>%
  filter(!is.na(log10_volume)) %>%
  group_by(genus) %>%
  summarise(
    size_logvol = median(log10_volume, na.rm = TRUE),
    phylum_heim = first(phylum), 
    class_heim = first(class)
  ) %>%
  ungroup()

dat <- dat %>%
  inner_join(heim_sizes, by = "genus")



#adding Heim et al. stages

dat$genus <- sub("_.*", "", dat$suggested_split_taxon)

genus_dat <- dat %>%
  group_by(genus) %>%
  summarise(
    FAD = max(max_ma, na.rm = TRUE),
    LAD = min(min_ma, na.rm = TRUE)
  )

stages <- stages_heim %>%
  rename(
    stage = Int,
    stage_name = Int_name,
    stage_old = old,
    stage_young = young
  ) %>%
  arrange(stage)   # Stage 1 = oldest, 89 = youngest

#stage_old >= FAD > stage_young

assign_FAD_stage <- function(fad_value) {
  match <- stages %>%
    filter(stage_old >= fad_value,
           fad_value > stage_young) %>%
    slice_max(stage)
  if(nrow(match) == 0) return(NA)
  return(match$stage)
}

#stage_old > LAD >= stage_young

assign_LAD_stage <- function(lad_value) {
  match <- stages %>%
    filter(stage_old > lad_value,
           lad_value >= stage_young) %>%
    slice_min(stage)  
  if(nrow(match) == 0) return(NA)
  return(match$stage)
}

genus_stg <- genus_dat %>%
  rowwise() %>%
  mutate(
    FAD_stage = assign_FAD_stage(FAD),
    LAD_stage = assign_LAD_stage(LAD)
  ) %>%
  ungroup()

last_stage <- max(genus_stg$LAD_stage, na.rm = TRUE)

dat <- dat %>%
  left_join(genus_stg, by = "genus")



#adding eco simple and small/large

dat <- dat %>%
  mutate(
    
    #motility
    motility_simplified = case_when(
      # Motile 
      str_detect(motility, "actively mobile|fast-moving|slow-moving|passively mobile|facultatively mobile") ~ "Motile",
      
      # Stationary
      str_detect(motility, "stationary|attached") ~ "Stationary",
      
      TRUE ~ "Other"
    ),
    
    #diet
    diet_simplified = case_when(
      # Carnivores (incl. microcarnivore, piscivore, durophage and parasites)
      str_detect(diet, "carnivore|piscivore|durophage|microcarnivore|parasite|parasitic") ~ "Carnivore",
      
      # Grazers / Herbivores / Browsers
      str_detect(diet, "grazer|herbivore|browser") ~ "Grazer",
      
      # Detritivores / Deposit feeders / Coprophages
      str_detect(diet, "detritivore|deposit feeder|coprophage") ~ "Detritivore",
      
      # Suspension feeders
      str_detect(diet, "suspension feeder") ~ "Suspension feeder",
      
      # Everything else (omnivores, chemo/photo-symbiotic)
      TRUE ~ "Other"
    ),
    
    #life habit
    life_habit_simplified = case_when(
      # Nektonic
      str_detect(life_habit, "nektonic|nektobenthic") ~ "Nektonic",
      
      # Infaunal
      str_detect(life_habit, "infaunal|deep infaunal|semi-infaunal|shallow infaunal") ~ "Infaunal",
      
      # Epifaunal
      str_detect(life_habit, "epifaunal|boring|upper-level epifaunal|low-level epifaunal|intermediate-level epifaunal") ~ "Epifaunal",
      
      TRUE ~ "Other"
    ))

# global median
global_median <- median(dat$size_logvol, na.rm = TRUE)

#small/large
dat <- dat %>%
  mutate(
    size_class = ifelse(size_logvol < global_median, "small", "large")
  )
#occurrence count 
occ_counts <- dat %>%
  +     count(genus, name = "occurrence_count")
dat <- dat %>%
  +     left_join(occ_counts, by = "genus")

saveRDS(dat, "dat.rds")


############################Step 2: data structure/sampling bias############

#scatterplot logsize/logocc

genus_scatter_data <- dat %>%
  group_by(genus) %>%  
  summarise(
    size_logvol = median(size_logvol, na.rm = TRUE),
    occurrences = n()
  ) %>%
  filter(!is.na(size_logvol))

mod_bias <- lm(log(occurrences) ~ size_logvol, data = genus_scatter_data)

summary(mod_bias)

ggplot(genus_scatter_data, aes(x = size_logvol, y = log(occurrences))) +
  geom_point(alpha = 0.2, color = "darkseagreen4") +
  geom_smooth(method = "lm", color = "grey33", linewidth = 0.6) +
  theme_minimal() +
  labs(
    title = "Body Size vs. Sampling Frequency",
    x = "Body Size (log-vol)",
    y = "Log(Occurrences)"
  )

#size distribution
plot_size_distribution <- function(data_input, factor_col, title_text, colors, top_n = NULL) {
  
  plot_data <- data_input %>%
    filter(!is.na(.data[[factor_col]])) %>%
    filter(.data[[factor_col]] != "Other")
  
  if (!is.null(top_n)) {
    top_categories <- plot_data %>%
      count(.data[[factor_col]], sort = TRUE) %>%
      slice_head(n = top_n) %>%
      pull(.data[[factor_col]])
    
    plot_data <- plot_data %>%
      filter(.data[[factor_col]] %in% top_categories)
  }
  
  global_median <- median(data_input$size_logvol, na.rm = TRUE)
  
  ggplot(plot_data, aes(x = reorder(.data[[factor_col]], size_logvol, FUN = median), 
                        y = size_logvol, 
                        fill = .data[[factor_col]])) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = colors, guide = "none") + 
    geom_hline(yintercept = global_median, color = "grey50", linetype = "dashed", linewidth = 1) +
    
    coord_cartesian(ylim = range(data_input$size_logvol, na.rm = TRUE)) +
    
    labs(
      title = title_text,
      y = "Body Size (log-vol)",
      x = "" 
    ) +
    theme_minimal(base_size = 12) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      panel.grid.major.x = element_blank(), 
      plot.title = element_text(face = "bold", size = 12) 
    )
  
  return(g)
}

plot_size_distribution(dat, "diet_simplified", 
                       "Diet",
                       c("lightsalmon4", "darkseagreen4", "palegreen3", "papayawhip"))

plot_size_distribution(dat, "life_habit_simplified", 
                       "Life Habit", 
                       c("lightsalmon4", "palegreen3", "papayawhip"))

plot_size_distribution(dat, "motility_simplified", 
                       "Motility", 
                       c("darkseagreen4", "papayawhip"))

#table for exact numbers
get_summary_stats <- function(data, group_var) {
  data %>%
    filter(!is.na(.data[[group_var]]), .data[[group_var]] != "Other") %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),                                      
      Median = median(size_logvol, na.rm = TRUE),
      Mean = mean(size_logvol, na.rm = TRUE),       
      SD = sd(size_logvol, na.rm = TRUE),           
      Min = min(size_logvol, na.rm = TRUE),
      Max = max(size_logvol, na.rm = TRUE)
    ) %>%
    arrange(desc(Median))
}

cat("GLOBAL MEDIAN:", median(dat$size_logvol, na.rm = TRUE), "\n\n")

print(get_summary_stats(dat, "diet_simplified"))

print(get_summary_stats(dat, "life_habit_simplified"))

print(get_summary_stats(dat, "motility_simplified"))


############################# Step 3: Extinction Pattern ##########

data(stages)
run_chisq <- function(count_matrix) {
  test <- chisq.test(count_matrix)
  print(test)
  
  rates <- count_matrix[,1] / rowSums(count_matrix)
  cat("--- Extinction Rates ---\n")
  print(rates)
  
  if(test$p.value < 0.05) cat(">> Result: Significance difference!\n")
  else cat(">> Result: No significant difference.\n")
}

prep_dd <- function(res) {
  if(is.null(res)) return(NULL)
  res$bin <- as.numeric(rownames(res))
  res$age <- stages$mid[res$bin]
  return(res)
}

col_small <- "tan"
col_large <- "darkseagreen4"

#global (Raw + Subsampled)

global_median <- median(dat$size_logvol, na.rm = TRUE)
dat_global <- dat %>%
  mutate(size_class = ifelse(size_logvol <= global_median, "Small", "Large"))

small_raw <- divDyn(dat_global %>% filter(size_class == "Small"), tax = "genus", bin = "stg") %>% prep_dd()
large_raw <- divDyn(dat_global %>% filter(size_class == "Large"), tax = "genus", bin = "stg") %>% prep_dd()

# 3. Subsampling 
small_sub <- subsample(dat_global %>% filter(size_class == "Small"), 
                       tax = "genus", bin = "stg", type = "sqs", q = 0.8, iter = 50) %>% prep_dd()
large_sub <- subsample(dat_global %>% filter(size_class == "Large"), 
                       tax = "genus", bin = "stg", type = "sqs", q = 0.8, iter = 50) %>% prep_dd()

#congruency testing, test for both small and large
comparison_data_large <- data.frame(
  Stage = as.numeric(rownames(large_raw)),
  Raw_Extl = large_raw$extPC,
  Sub_Extl = large_sub$extPC
)
comparison_data_small <- data.frame(
  Stage = as.numeric(rownames(large_raw)),
  Raw_Exts = small_raw$extPC,
  Sub_Exts = small_sub$extPC
)

comparison_data_large <- na.omit(comparison_data_large)
comparison_data_small <- na.omit(comparison_data_small)

test_result_large <- cor.test(comparison_data_large$Raw_Extl, 
                              comparison_data_large$Sub_Extl, 
                              method = "spearman")

test_result_small <- cor.test(comparison_data_small$Raw_Exts, 
                              comparison_data_small$Sub_Exts, 
                              method = "spearman")
print(test_result_large)
print(test_result_small)

#plot
tsplot(stages, boxes = "sys", shading = "sys", xlim = c(520, 24), ylim = c(0, 1.0), 
       xlab = "Age (Ma)", ylab = "Extinction Rate (extPC)")
title(main = "Global Extinction: Small vs. Large", cex.main = 1.25)

lines(small_raw$age, small_raw$extPC, col = col_small, lty = 2, lwd = 2)
lines(large_raw$age, large_raw$extPC, col = col_large, lty = 2, lwd = 2)

lines(small_sub$age, small_sub$extPC, col = col_small, lty = 1, lwd = 3)
lines(large_sub$age, large_sub$extPC, col = col_large, lty = 1, lwd = 3)

legend("topright", legend = c("Small (Sub)", "Large (Sub)", "Small (Raw)", "Large (Raw)"),
       col = c(col_small, col_large, col_small, col_large),
       lty = c(1, 1, 2, 2), lwd = c(4, 4, 2, 2), bg = "white")

n_ext_s <- sum(small_raw$tExt, na.rm=T)
n_surv_s <- sum(small_raw$divBC, na.rm=T) 
n_ext_l <- sum(large_raw$tExt, na.rm=T)
n_surv_l <- sum(large_raw$divBC, na.rm=T)

mat_global <- matrix(c(n_ext_s, n_surv_s, n_ext_l, n_surv_l), nrow = 2, byrow = T,
                     dimnames = list(c("Small", "Large"), c("Extinct", "Survived")))

run_chisq(mat_global)

#eco triats

eco_traits <- c("diet_simplified", "motility_simplified", "life_habit_simplified") 

analyze_and_plot_trait <- function(df, trait_col, plot_title, my_colors) {
  
  dat_eco <- df %>% 
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "Other")
  categories <- sort(unique(dat_eco[[trait_col]]))
  
  tsplot(stages, boxes = "sys", shading = "sys", xlim = c(520, 24), ylim = c(0, 1.0),
         ylab = "Extinction Rate (extPC)")
  
  title(main = plot_title, cex.main = 1.25) 
  
  chi_names <- c()
  chi_ext <- c()
  chi_surv <- c()
  for(i in seq_along(categories)) {
    cat_name <- categories[i]
    color_for_cat <- my_colors[i]
    sub_dat <- dat_eco %>% filter(.data[[trait_col]] == cat_name)
    if(nrow(sub_dat) > 50) {
      res_eco <- divDyn(sub_dat, tax = "genus", bin = "stg")
      res_eco$age <- stages$mid[as.numeric(rownames(res_eco))]
      lines(res_eco$age, res_eco$extPC, col = color_for_cat, lwd = 3)
      chi_names <- c(chi_names, cat_name)
      col_ext <- if("tExt" %in% names(res_eco)) "tExt" else "ext"
      
      chi_ext <- c(chi_ext, sum(res_eco[[col_ext]], na.rm=T))
      chi_surv <- c(chi_surv, sum(res_eco$divBC, na.rm=T))
    } 
  }
  
  legend("topright", legend = categories, col = my_colors, lwd = 3, bg = "white", cex=0.8)
  
  if(length(chi_names) > 1) {
    mat_eco <- matrix(c(chi_ext, chi_surv), ncol = 2, 
                      dimnames = list(chi_names, c("Extinct", "Survived")))
    run_chisq(mat_eco)
  }
}

#diet
analyze_and_plot_trait(
  df = dat, 
  trait_col = "diet_simplified", 
  plot_title = "Extinction Rates by Diet", 
  my_colors = c("lightsalmon4", "darkseagreen4", "palegreen3", "papayawhip") 
)

#motility
analyze_and_plot_trait(
  df = dat, 
  trait_col = "motility_simplified", 
  plot_title = "Extinction Rates by Motility", 
  my_colors = c("darkseagreen4", "papayawhip") 
)

#life habit
analyze_and_plot_trait(
  df = dat, 
  trait_col = "life_habit_simplified", 
  plot_title = "Extinction Rates by Life Habit", 
  my_colors = c("lightsalmon4", "palegreen3", "papayawhip")
)


#GLM Extinction in Stage ~ body size

stgs <- function() {
  data(stages)
  
  stages_dat <- stages %>%
    rename(
      stage_id = stg,
      bottom = bottom,
      top = top,
      mid_age = mid
    ) %>%
    dplyr::select(stage_id, bottom, top, mid_age, sys)
  
  return(stages_dat)
}

longdat <- function(dat) {
  dat_clean <- dat %>%
    filter(is.finite(FAD_stage), is.finite(LAD_stage))
  
  dat_long_raw <- dat_clean %>%
    dplyr::select(phylum, genus, FAD_stage, LAD_stage, 
                  size_logvol, diet_simplified, motility_simplified, life_habit_simplified) %>% # occ_bin entfernt!
    filter(LAD_stage >= FAD_stage) %>% 
    rowwise() %>%
    mutate(stage_id = list(seq(FAD_stage, LAD_stage))) %>%
    unnest(stage_id) %>%
    ungroup()
  
  dat_genus <- dat_long_raw %>%
    group_by(stage_id, genus) %>%
    summarise(
      occ_raw = n(), # Zählt die Vorkommen pro Gattung und Stufe
      size_logvol = mean(size_logvol, na.rm = TRUE), 
      diet_simplified = first(diet_simplified),
      motility_simplified = first(motility_simplified),
      life_habit_simplified = first(life_habit_simplified),
      phylum = first(phylum),
      LAD_stage = first(LAD_stage), # occ_bin auch hier entfernt!
      .groups = "drop"
    ) %>%
    mutate(
      extinct_in_stage = if_else(stage_id == LAD_stage, 1L, 0L),
      occ_log = log1p(occ_raw) # Logarithmierte Occurrences pro Stufe. Das nutzen wir!
    )
  
  return(dat_genus)
}


run_stagewise_glm <- function(data, formula, target_var = "size_logvol") {
  
  results <- data %>%
    group_by(stage_id) %>%
    filter(n() > 10, n_distinct(extinct_in_stage) == 2) %>%
    summarise(
      model = list(glm(formula, family = binomial(), data = pick(everything()))),
      coef = tryCatch(coef(model[[1]])[target_var], error = function(e) NA),
      se   = tryCatch(summary(model[[1]])$coefficients[target_var, 2], error = function(e) NA),
      pval = tryCatch(summary(model[[1]])$coefficients[target_var, 4], error = function(e) NA),
      
      n_total = n(),
      n_ext = sum(extinct_in_stage == 1),
      .groups = "drop"
    ) %>%
    filter(!is.na(coef))
  
  return(results)
}

summ_results <- function(results, model_name = "Model") {
  
  n_total_stages <- nrow(results)
  sig_results <- results %>% filter(pval < 0.05)
  n_sig <- nrow(sig_results)
  
  cat("\n--- Summary for:", model_name, "---\n")
  cat("Total stages analyzed:", n_total_stages, "\n")
  cat("Significant stages (p < 0.05):", n_sig, paste0("(", round((n_sig / n_total_stages) * 100, 1), "%)\n"))
  
  if(n_sig > 0) {
    n_sig_neg <- sum(sig_results$coef < 0) 
    n_sig_pos <- sum(sig_results$coef > 0) 
    
    cat(" -> Negative effect (Variable protects):", n_sig_neg, "\n")
    cat(" -> Positive effect (Variable kills):", n_sig_pos, "\n")
  }
  cat("----------------------------------\n")
  
  return(invisible(sig_results))
}

# plot
plot_stagewise_effect <- function(results, stages_ref, 
                                  main_title = "", 
                                  ylab = "Log-Odds (Size Effect)",
                                  col = "black",
                                  ylim = c(-1, 1)){
  
  plot_dat <- results %>%
    left_join(stages_ref, by = "stage_id")
  max_val <- max(abs(c(plot_dat$coef + 1.96*plot_dat$se, plot_dat$coef - 1.96*plot_dat$se)), na.rm = TRUE)
  limit <- max(1.5, max_val * 1.1) 
  tsplot(stages, boxes = "sys", shading = "sys", 
         xlim = c(520, 24),  
         ylim = c(-1, 1),
         xlab = "Age (Ma)",
         ylab = ylab,
         prop = 0.05)
  
  title(main = main_title, line = 1.5, cex.main = 1.2)
  abline(h = 0, lty = 2, col = "grey50", lwd = 1.5)
  
  segments(x0 = plot_dat$mid_age, 
           y0 = plot_dat$coef - 1.96 * plot_dat$se,
           x1 = plot_dat$mid_age, 
           y1 = plot_dat$coef + 1.96 * plot_dat$se,
           col = "grey30", lwd = 1.2)
  
  lines(plot_dat$mid_age, plot_dat$coef, type = "l", lwd = 2.5, col = col)
  
  points(plot_dat$mid_age, plot_dat$coef, pch = 21, bg = "white", col = col, cex = 1.0)
  sig_points <- plot_dat %>% filter(pval < 0.05)
  if(nrow(sig_points) > 0) {
    points(sig_points$mid_age, sig_points$coef, pch = 21, bg = col, col = "black", cex = 1.3)
  }
  
  legend("topright", legend = c("Significant", "Not significant"), 
         pch = 21, pt.bg = c(col, "white"), col = c("black", col), bty = "n", cex = 0.8)
}

stages_dat <- stgs()

dat_long<- longdat (dat) 
res_size <- run_stagewise_glm(dat_long, extinct_in_stage ~ size_logvol, target_var = "size_logvol")

plot_stagewise_effect(res_size, stages_dat, 
                      main_title = "Extinction ~ log(Body Size)",
                      col = "darkseagreen4")
summ_results(res_size, "GLM bodysize")



###########################Step 4: Duration Patterns#####################


#forest plot/multiple model
library(MASS)
library(broom)

model_data <- dat %>%
  group_by(genus) %>%
  summarise(
    duration_ma = max(max_ma, na.rm = TRUE) - min(min_ma, na.rm = TRUE),
    size_logvol = median(size_logvol, na.rm = TRUE),
    occurrences = n(), 
    diet = first(diet_simplified),
    motility = first(motility_simplified),
    life_habit = first(life_habit_simplified)
  ) %>%
  filter(duration_ma > 0, 
         !is.na(size_logvol), 
         !is.na(diet), diet != "Other", 
         !is.na(motility), motility != "Other",
         !is.na(life_habit), life_habit !="Other")

full_model_duration <- lm(duration_ma ~ size_logvol + log(occurrences) + diet + motility + life_habit, 
                          data = model_data)
#AIC result
best_model_duration <- stepAIC(full_model_duration, direction = "both", trace = 0)

#result
summary(best_model_duration)

# plots
model_results <- tidy(best_model_duration, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term_clean = case_when(
    term == "size_logvol" ~ "Body Size (log-vol)",
    term == "log(occurrences)" ~ "Sampling Intensity",
    term == "dietDetritivore" ~ "Detritivore",
    term == "dietGrazer" ~ "Grazer",
    term == "dietSuspension feeder" ~ "Suspension Feeder",
    term == "motilityStationary" ~ "Stationary",
    term == "life_habitInfaunal" ~ "Infaunal",
    term == "life_habitNektonic" ~ "Nektonic",
    term == "life_habitOther" ~ "Other",
    TRUE ~ term
  )) %>%
  
  #sort
  mutate(term_ordered = factor(term_clean, levels = rev(term_clean))) %>%
  mutate(Effect = ifelse(estimate > 0, "Protective (Longer Life)", "Risk (Shorter Life)"))

model_results_grouped <- model_results %>%
  mutate(Group = case_when(
    term %in% c("size_logvol", "log(occurrences)") ~ "Controls (Continuous)",
    grepl("diet", term) ~ "Diet (ref: Carnivore)",
    grepl("motility", term) ~ "Motility (ref: Motile)",
    grepl("life_habit", term) ~ "Life Habit (ref: Epifaunal)"
  ))

ggplot(model_results_grouped, aes(x = estimate, y = term_clean, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3.5) + 
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, size = 0.8) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(
    title = "Drivers of Evolutionary Duration",
    subtitle = "Multiple Linear Regression Model (AIC-selected)",
    x = "Effect on Duration (Million Years)",
    y = ""
  ) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 10, hjust = 0),
    strip.background = element_rect(fill = "gray95")
  ) +
  scale_color_manual(values = c("Protective (Longer Life)" = "darkseagreen4", "Risk (Shorter Life)" = "lightsalmon4"))

#partial corr
library(car)

mod_occ <- lm(duration_ma ~ log(occurrences), data = model_data)
print(summary(mod_occ))

mod_size <- lm(duration_ma ~ size_logvol, data = model_data)
print(summary(mod_size))

mod_both <- lm(duration_ma ~ log(occurrences) + size_logvol, data = model_data)
print(summary(mod_both))

avPlot(mod_both, 
       variable = "size_logvol", 
       pch = 16,                 
       col = rgb(0.5, 0.7, 0.5, 0.5), 
       col.lines = "grey45",        
       lwd = 2,    
       id = FALSE,
       main = "Partial Correlation: Body Size vs. Duration\n(Controlling for Sampling Intensity)",
       xlab = "Body Size (log-vol) | given log(Occurrences)",
       ylab = "Evolutionary Duration (Ma) | given log(Occurrences)")