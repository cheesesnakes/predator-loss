pacman::p_load(dplyr, tidyr, here, stringr, purrr, flextable)

# import data

sites <- read.csv(here("data/chapter-3_sites.csv"))
plots <- read.csv(here("data/chapter-3_plots.csv"))
rugosity <- read.csv(here("data/chapter-3_rugosity.csv"))
benthic <- read.csv(here("data/chapter-3_benthic.csv"))


samples <- read.csv(here("data/chapter-3_samples.csv"))
fish <- read.csv(here("data/chapter-3_individuals.csv"))
groups <- read.csv(here("data/chapter-3_groups.csv"))

behaviors <- read.csv(here("data/chapter-3_behaviors.csv"))
boris <- read.csv(here("data/chapter-3_behaviors_boris.tsv"), sep = "\t", header = TRUE)

# Check and create directories

dir.create(here("outputs"), showWarnings = FALSE)
dir.create(here("outputs/errors"), showWarnings = FALSE)

# PREDICTORS -------------------------------------------------------------

# site metadata

sites <- sites %>%
  select(deployment.id, location, site, protection)

# plot metadata

resident_damsels <- plots %>%
  select(resident_damsels) %>%
  distinct() %>%
  pull()

damsels <- plots %>%
  select(plot_id, damsel_n) %>%
  # damsel presence
  mutate(damsel = ifelse(damsel_n > 0, 1, 0)) %>%
  select(-damsel_n) %>%
  distinct() %>%
  filter(!is.na(damsel))

# benthic cover

category_names <- read.csv(here("data/benthic_categories.csv"))

benthic <- benthic %>%
  left_join(category_names, by = "category") %>%
  select(plot_id, names, count) %>%
  mutate(prop = count / 200) %>%
  select(plot_id, names, prop) %>%
  filter(prop > 0) %>%
  group_by(plot_id, names) %>%
  summarise(cover = sum(prop)) %>%
  ungroup() %>%
  pivot_wider(names_from = names, values_from = cover) %>%
  replace(is.na(.), 0)

# rugosity data

rugosity <- rugosity %>%
  # create plot_id join with deployment.id
  rename(deployment.id = plot_id) %>%
  mutate(plot_id = paste0(deployment.id, "", plot.id)) %>%
  # calculate rugosity
  mutate(rugosity = 1 - (d_reef_cm / chain.length.cm))

rugosity <- rugosity %>%
  select(plot_id, transect.no, rugosity) %>%
  pivot_wider(names_from = transect.no, values_from = rugosity, names_prefix = "sample_") %>%
  replace(is.na(.), 0)


# plot-level rugosity and variation

rugosity <- rugosity %>%
  mutate(
    mean = rowMeans(select(., starts_with("sample_")), na.rm = TRUE),
    se = apply(select(., starts_with("sample_")), 1, sd) / sqrt(3)
  )

moments_rugosity <- c(mean(rugosity$mean, na.rm = TRUE), sd(rugosity$mean, na.rm = TRUE))

rugosity <- rugosity %>%
  mutate(
    mean = scale(mean)[, 1],
    se = se / moments_rugosity[2]
  ) %>%
  select(plot_id, mean, se) %>%
  rename(rugosity_mean = mean, rugosity_se = se)

# predator presence

predator <- samples %>%
  mutate(plot_id = paste0(deployment_id, "", plot_id)) %>%
  mutate(sample_id = paste0(plot_id, "", sample_id)) %>%
  # code predator presence
  mutate(predator = ifelse(predator.present == "yes", 1, 0)) %>%
  select(plot_id, predator) %>%
  filter(predator == 1) %>%
  distinct()

# make predictors

predictors <- data.frame(plot_id = unique(rugosity$plot_id)) %>%
  mutate(deployment.id = str_remove(plot_id, "[0-9]")) %>%
  left_join(sites, by = "deployment.id") %>%
  left_join(benthic, by = "plot_id") %>%
  left_join(rugosity, by = "plot_id") %>%
  left_join(predator, by = "plot_id") %>%
  arrange(deployment.id, plot_id) %>%
  filter(deployment.id != "FBWT-P") %>%
  group_by(deployment.id) %>%
  # impute missing values with mean
  mutate_at(c("biomass", "substrate", "coral", "sponge"), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
  ungroup() %>%
  left_join(damsels, by = "plot_id") %>%
  replace(is.na(.), 0) %>%
  select(-deployment.id)

# Sampled Individuals -------------------------------------------------------------

fish <- fish %>%
  mutate(group = ifelse(group_id != "", 1, 0)) %>%
  mutate(size = gsub(" ", "", size))

# clean size

fish <- fish %>%
  # trim spaces
  mutate(size = gsub(" ", "", size)) %>%
  # replace 32 with 30
  mutate(size = gsub("32", "30", size)) %>%
  # make 0 - 10 into <10
  mutate(size = gsub("0-10", "<10", size)) %>%
  # make 10–20 into 10-20
  mutate(size = gsub("10–20", "10-20", size)) %>%
  # make 20–30 into 20-30
  mutate(size = gsub("20–30", "20-30", size)) %>%
  # arrange as factor <10, 10-20, 20-30, 30-40, 40-50
  mutate(size = factor(size, levels = c("<10", "10-20", "20-30", "30-40", "40-50")))

# remove OTU and unidentified species

fish <- fish %>%
  filter(!grepl("OTU", spp) & !grepl("spp", spp) & spp != "") %>%
  # remove green sea turtle
  filter(!grepl("Chelonia mydas", spp)) %>%
  # remove leading and trailing spaces
  mutate(spp = str_trim(spp))


fish <- fish %>%
  select(c(ind_id, sample_id, spp, group, size))

samples <- samples %>%
  select(c(sample_id, plot_id, deployment_id)) %>%
  mutate(
    plot_id = paste0(deployment_id, "", plot_id),
    sample_id = paste0(plot_id, "", sample_id)
  ) %>%
  select(-deployment_id)

fish <- fish %>%
  left_join(samples, by = "sample_id") %>%
  select(-c(sample_id))

# Make plot-level abudance data

plot_abundance <- fish %>%
  group_by(plot_id, spp) %>%
  summarise(abundance = n()) %>%
  ungroup() %>%
  mutate(resident = ifelse(spp %in% resident_damsels, 1, 0))

plot_size <- fish %>%
  group_by(plot_id, size) %>%
  summarise(abundance = n())


# attach foraging guild information

traits <- read.csv(here("data", "fish-traits.csv"))

phylo <- traits %>%
  select(c(Order, Family, Genus)) %>%
  filter(Family != "") %>%
  distinct()

traits <- traits %>%
  select(c(Family, Genus, Species, Length..cm., Height, Feeding.Guild, Habitat, Water.Column)) %>%
  mutate(species = paste(Genus, Species)) %>%
  rename(
    length = Length..cm.,
    height = Height, guild = Feeding.Guild, habitat = Habitat, column = Water.Column
  ) %>%
  select(-c(Species)) %>%
  mutate_if(~ is.character(.), ~ str_trim(.)) %>%
  mutate(Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus)) %>%
  mutate(guild = str_to_sentence(guild))


# Helper function to get mode
get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) {
    return(NA_character_)
  }
  tbl <- table(x)
  max_freq <- max(tbl)
  mode_candidates <- names(tbl)[tbl == max_freq]
  return(mode_candidates[1]) # return the first if there's a tie # nolint: return_linter.
}

# Trait imputation function
trait_mode <- function(family = NA, genus = NA, trait = NA) {
  # Try genus-level
  trait_vals <- traits %>%
    filter(Genus == genus) %>%
    pull(!!sym(trait))

  mode_val <- get_mode(trait_vals)

  # If genus-level fails, try family-level
  if (is.na(mode_val)) {
    trait_vals <- traits %>%
      filter(Family == family) %>%
      pull(!!sym(trait))
    mode_val <- get_mode(trait_vals)
  }

  return(mode_val)
}

# Main species data processing

fish_traits <- data.frame(species = unique(fish$spp)) %>%
  left_join(traits, by = "species") %>%
  mutate(
    Genus = ifelse(is.na(Genus), str_split(species, " ", simplify = TRUE)[, 1], Genus)
  ) %>%
  left_join(phylo, by = "Genus", relationship = "many-to-many") %>%
  mutate(
    Family = ifelse(is.na(Family.x), Family.y, Family.x),
    guild = ifelse(
      is.na(guild) | guild == "Unknown",
      map2_chr(Family, Genus, ~ trait_mode(family = .x, genus = .y, trait = "guild")),
      guild
    ),
    guild = na_if(guild, "NULL") # clean up any residual "NULL"
  ) %>%
  mutate(Family = ifelse(species == "Zebrasoma rostratum", "Acanthuridae", Family))

# Get family information

families <- fish_traits %>%
  select(Family, species) %>%
  distinct()

# add to fish data

fish <- fish %>%
  left_join(families, by = c("spp" = "species"))

colnames(fish) <- str_to_lower(colnames(fish))

# Mark residents

fish <- fish %>%
  mutate(resident = ifelse(spp %in% resident_damsels, 1, 0))

# calculate behavioral traits

# make time series data

responses <- behaviors %>%
  # make all columns numeric except ind_id
  group_by(ind_id) %>%
  mutate_if(is.character, as.numeric) %>%
  ungroup() %>%
  # fill NAS
  replace(is.na(.), 0) %>%
  # calculate foraging time
  mutate(foraging_time = forage_end - forage_start) %>%
  # calculate vigilance time
  mutate(vigilance_time = vigil_end - vigil_start) %>%
  # calculate move time
  mutate(move_time = move_end - move_start) %>%
  select(ind_id, obs_id, foraging_time, vigilance_time, move_time, bites, escape)

aggression <- behaviors %>%
  group_by(ind_id) %>%
  pivot_longer(cols = c(con_ag, het_ag, pred_ag), names_to = "modifier", values_to = "aggression") %>%
  select(ind_id, obs_id, aggression, modifier) %>%
  drop_na(aggression)

responses <- responses %>%
  left_join(aggression) %>%
  mutate(aggression = ifelse(is.na(aggression), 0, aggression))

# clean boris data

boris <- boris %>%
  select(Observation.id, Coding.duration, Behavior, Modifier..1, Behavior.type, Duration..s., Start..s., Stop..s.) %>%
  # make all colnames lower case
  rename_all(tolower) %>%
  # remove '.'
  rename_all(~ gsub("\\.", "", .)) %>%
  group_by(observationid, behavior) %>%
  mutate(count = n()) %>%
  distinct(observationid, codingduration, behavior, behaviortype, modifier1, starts, stops, durations, count)

boris_state <- boris %>%
  filter(behaviortype == "STATE") %>%
  select(observationid, codingduration, behavior, durations) %>%
  group_by(observationid, codingduration, behavior) %>%
  summarise(durations = sum(durations)) %>%
  ungroup() %>%
  pivot_wider(names_from = behavior, values_from = durations) %>%
  replace(is.na(.), 0) %>%
  rename(
    ind_id = observationid,
    foraging_time = Foraging,
    vigilance_time = Vigilant,
    move_time = Moving
  ) %>%
  # change types
  mutate(
    foraging_time = as.numeric(foraging_time),
    vigilance_time = as.numeric(vigilance_time),
    move_time = as.numeric(move_time)
  )

boris_point <- boris %>%
  filter(behaviortype == "POINT") %>%
  group_by(observationid, codingduration, behavior, modifier1) %>%
  summarise(count = last(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = behavior, values_from = count) %>%
  replace(is.na(.), 0) %>%
  rename(
    ind_id = observationid,
    bites = Bite,
    escape = Escape,
    aggression = Agression,
    modifier = modifier1
  ) %>%
  # change types
  mutate(
    bites = as.numeric(bites),
    escape = as.numeric(escape),
    aggression = as.numeric(aggression)
  ) %>%
  group_by(ind_id) %>%
  summarise(
    bites = sum(bites),
    escape = sum(escape),
    aggression = sum(aggression),
    modifier = last(modifier),
    codingduration = max(codingduration)
  ) %>%
  ungroup()

responses_boris <- boris_state %>%
  left_join(boris_point)

# join responses

responses <- responses %>%
  full_join(responses_boris) %>%
  # fill all numeric columns with 0
  mutate_if(is.numeric, ~ replace_na(., 0))

responses <- responses %>%
  select(-c(obs_id, modifier)) %>%
  group_by(ind_id) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  ungroup()

# Fix coding duration erros from boris

responses <- responses %>%
  mutate(codingduration = foraging_time + vigilance_time + move_time) %>%
  # fix coding durations longer than 120 s (bug in boris)
  mutate(
    foraging = ifelse(codingduration > 150, (foraging_time * 120 / codingduration), foraging_time),
    vigilance = ifelse(codingduration > 150, (vigilance_time * 120 / codingduration), vigilance_time),
    move = ifelse(codingduration > 150, (move_time * 120 / codingduration), move_time),
    codingduration = ifelse(codingduration > 150, 120, codingduration)
  )


## check for errors

fish %>%
  filter(!ind_id %in% responses$ind_id) %>%
  write.csv(here("outputs/errors/individuals_missing_behaviour.csv"), row.names = FALSE)

responses %>%
  filter(!ind_id %in% fish$ind_id) %>%
  write.csv(here("outputs/errors/behaviour_missing_individuals.csv"), row.names = FALSE)

fish <- fish %>%
  # remove individuals with no responses
  filter(ind_id %in% responses$ind_id)

responses <- responses %>%
  # remove individuals with no fish
  filter(ind_id %in% fish$ind_id) %>%
  filter(!is.na(foraging) & !is.na(vigilance) & !is.na(move))

## make instantaeneous scans

scans <- behaviors %>%
  select(ind_id, obs_id, forage_start, forage_end, vigil_start, vigil_end, move_start, move_end) %>%
  group_by(ind_id, obs_id) %>%
  pivot_longer(cols = c(forage_start, forage_end, vigil_start, vigil_end, move_start, move_end), names_to = "behavior", values_to = "time") %>%
  filter(!is.na(time)) %>%
  mutate(
    pos = str_split(behavior, "_", simplify = TRUE)[, 2],
    behavior = str_split(behavior, "_", simplify = TRUE)[, 1]
  ) %>%
  pivot_wider(names_from = pos, values_from = time) %>%
  mutate(scan = ifelse(start <= 30 & !(end <= 10), 1, 0)) %>%
  filter(scan == 1) %>%
  ungroup() %>%
  select(ind_id, behavior, scan) %>%
  pivot_wider(names_from = behavior, values_from = scan, values_fill = list(scan = 0))

scans_boris <- boris %>%
  filter(behaviortype == "STATE") %>%
  select(observationid, behavior, starts, stops) %>%
  group_by(observationid) %>%
  mutate(duration = stops - starts) %>%
  arrange(observationid, starts) %>%
  # make duration cumulative
  mutate(cum_duration = cumsum(duration)) %>%
  mutate(starts = ifelse(cum_duration == min(cum_duration), 0, cum_duration - duration),
         stops = cum_duration,
         scan = ifelse(starts <= 30 & !(stops <= 10), 1, 0)) %>%
  filter(scan == 1) %>%
  select(observationid, behavior, duration, scan) %>%
  mutate(scan = ifelse(duration < max(duration), 0, scan)) %>%
         filter(scan == 1) %>%
         select(-c(duration)) %>%
         pivot_wider(names_from = behavior, values_from = scan, values_fill = list(scan = 0))

colnames(scans) <- c("ind_id", "vigilance", "foraging", "moving")
colnames(scans_boris) <- c("ind_id", "foraging", "moving", "vigilance")

scans <- rbind(scans, scans_boris)

# Guild membership

guild <- fish_traits %>%
  select(species, guild) %>%
  mutate(
    FG1 = str_split(guild, ", ", simplify = T)[, 1],
    FG2 = str_split(guild, ", ", simplify = T)[, 2]
  ) %>%
  select(-c(guild)) %>%
  pivot_longer(cols = c(FG1, FG2), names_to = "type", values_to = "guild") %>%
  select(-c(type)) %>%
  mutate(guild = str_to_sentence(str_trim(guild))) %>%
  distinct() %>%
  filter(guild != "") %>%
  mutate(value = 1) %>%
  complete(species, guild, fill = list(value = 0)) %>%
  group_by(species) %>%
  reframe(
    value = value / sum(value),
    guild = guild
  ) %>%
  pivot_wider(names_from = guild, values_from = value)

# multiple guilds

multiguild <- guild[which(rowSums(guild[, -1] > 0) > 1), ] %>%
  select(species) %>%
  filter(species %in% fish$spp) %>%
  rename(Species = species)

print(paste0("Number of species with multiple guilds: ", nrow(multiguild)))

multiguild <- multiguild %>%
  flextable() %>%
  autofit()

save_as_image(
  multiguild,
  path = here("figures", "multiguild.png"),
  type = "png",
  width = 12,
  height = 8
)

# make multiple guild into single

guild <- guild %>%
  mutate(
    Herbivore = ifelse(Herbivore > 0, 1, 0),
    Piscivore = ifelse(Piscivore > 0, 1, 0),
    Invertivore = ifelse(Invertivore == 1, 1, 0),
  ) %>%
  pivot_longer(cols = c(Herbivore, Piscivore, Invertivore), names_to = "guild", values_to = "value") %>%
  filter(value > 0) %>%
  select(-c(value))

# remove temporary objects

rm(
 benthic, boris_point, boris, boris_state, category_names, damsels,
 predator, sites, aggression, behaviors, groups, plots, responses_boris, samples, scans_boris
)

# final data

head(responses)
head(predictors)
head(fish)
head(traits)
head(guild)
