pacman::p_load(dplyr, tidyr, here, stringr, purrr)

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

# PREDICTORS -------------------------------------------------------------

# site metadata

sites <- sites%>%
    select(deployment.id, location, site, protection)

# plot metadata

damsels <- plots%>%
    select(plot_id, damsel_n)%>%
    # damsel presence
    mutate(damsel = ifelse(damsel_n > 0, 1, 0))%>%
    select(-damsel_n)%>%
    distinct()

   
# benthic cover 

category_names = read.csv(here("data/benthic_categories.csv"))

benthic <- benthic%>%
    left_join(category_names, by = "category")%>%
    select(plot_id, names, count)%>%
    mutate(prop = count/200)%>%
    select(plot_id, names, prop)%>%
    filter(prop>0)%>%
    group_by(plot_id, names)%>%
    summarise(cover = sum(prop))%>%
    ungroup()%>%
    pivot_wider(names_from = names, values_from = cover)%>%
    replace(is.na(.), 0)

# rugosity data

rugosity <- rugosity%>%
    # create plot_id join with deployment.id
    rename(deployment.id = plot_id)%>%
    mutate(plot_id = paste0(deployment.id, "", plot.id))%>%
    # calculate rugosity
    mutate(rugosity = 1 - (d_reef_cm/chain.length.cm))
    
mean_rugosity <- rugosity%>%
    group_by(plot_id)%>%
    summarise(rugosity = mean(rugosity, na.rm = TRUE))%>%
    ungroup()

rugosity <- rugosity%>%
    select(plot_id, transect.no, rugosity)%>%
    pivot_wider(names_from = transect.no, values_from = rugosity, names_prefix = "sample_")%>%
    replace(is.na(.), 0)

# predator presence

predator <- samples%>%
mutate(plot_id = paste0(deployment_id, "", plot_id))%>%
mutate(sample_id = paste0(plot_id, "", sample_id))%>%
# code predator presence
mutate(predator = ifelse(predator.present == "yes", 1, 0))%>%
select(plot_id, predator)%>%
filter(predator == 1)%>%
distinct()

# make predictors

predictors <- data.frame(plot_id  = unique(rugosity$plot_id))%>%
mutate(deployment.id = str_remove(plot_id, "[0-9]"))%>%
left_join(sites, by = "deployment.id")%>%
left_join(benthic, by = "plot_id")%>%
left_join(mean_rugosity, by = "plot_id")%>%
left_join(predator, by = "plot_id")%>%
filter(!is.na(rugosity))%>%
arrange(deployment.id, plot_id)%>%
filter(deployment.id != "FBWT-P")%>%
group_by(deployment.id)%>%
# impute missing values with mean
mutate_at(c("resource", "substrate", "coral", "sponge"), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))%>%
ungroup()%>%
left_join(damsels, by = "plot_id")%>%
replace(is.na(.), 0)%>%
select(-deployment.id)

# Sampled Individuals -------------------------------------------------------------

fish <- fish%>%
mutate(group = ifelse(group_id != "", 1, 0))%>%
mutate(size = gsub(" ", "", size))

# clean size  

fish <- fish%>%
  # trim spaces
  mutate(size = gsub(" ", "", size))%>%
  # replace 32 with 30
  mutate(size = gsub("32", "30", size))%>%
  # make 0 - 10 into <10
  mutate(size = gsub("0-10", "<10", size))%>%
  # make 10–20 into 10-20
  mutate(size = gsub("10–20", "10-20", size))%>%
  # make 20–30 into 20-30
  mutate(size = gsub("20–30", "20-30", size))%>%
  # arrange as factor <10, 10-20, 20-30, 30-40, 40-50
  mutate(size = factor(size, levels = c("<10", "10-20", "20-30", "30-40", "40-50")))

# remove OTU and unidentified species

fish <- fish%>%
filter(!grepl("OTU", spp) & !grepl("spp", spp) & spp != "")%>%
# remove green sea turtle
filter(!grepl("Chelonia mydas", spp))%>%
#remove leading and trailing spaces
mutate(spp = str_trim(spp))


fish <- fish%>%
select(c(ind_id, sample_id, spp, group, size))

samples <- samples%>%
select(c(sample_id, plot_id, deployment_id))%>%
mutate(plot_id = paste0(deployment_id, "", plot_id),
sample_id = paste0(plot_id, "", sample_id))%>%
select(-deployment_id)

fish <- fish%>%
left_join(samples, by = "sample_id")%>%
select(-c(sample_id))

colnames(fish) <- str_to_lower(colnames(fish))

# attach foraging guild information

traits <- read.csv("data/transect-fish-traits.csv")

traits <- traits%>%
select(c(Genus, Species, Feeding.Guild))%>%
# join genus and species
mutate(spp = paste0(Genus, " ", Species))%>%
# split feeding guild into FG1 and FG2
mutate(Feeding.Guild = gsub(" ", "", Feeding.Guild))%>%
separate(Feeding.Guild, into = c("FG1", "FG2"), sep = ",")%>%
pivot_longer(cols = c(FG1, FG2), names_to = "FG", values_to = "Feeding.Guild")%>%
select(spp, Feeding.Guild)%>%
filter(Feeding.Guild != "Unknown", !is.na(Feeding.Guild))%>%
mutate(Feeding.Guild = str_to_sentence(Feeding.Guild),
Feeding.Guild =  ifelse(Feeding.Guild == "Corallivore", "Coralivore", Feeding.Guild),
Feeding.Guild =  ifelse(Feeding.Guild == "Planktivores", "Planktivore", Feeding.Guild))%>%
mutate(value = 1)%>%
complete(spp, Feeding.Guild, fill = list(value = 0))%>%
group_by(spp)%>%
mutate(value =  value/sum(value))%>%
pivot_wider(names_from = "Feeding.Guild", values_from = "value", values_fn = last)

# calculate behavioral traits

responses <- behaviors%>%
# make all columns numeric except ind_id
group_by(ind_id)%>%
mutate_if(is.character, as.numeric)%>%
ungroup()%>%
#fill NAS
replace(is.na(.), 0)%>%
# calculate foraging time
mutate(foraging_time = forage_end - forage_start)%>%
# calculate vigilance time
mutate(vigilance_time = vigil_end - vigil_start)%>%
# calculate move time
mutate(move_time = move_end - move_start)%>%
select(ind_id, obs_id, foraging_time, vigilance_time, move_time, bites, escape)

aggression <- behaviors%>%
group_by(ind_id)%>%
pivot_longer(cols = c(con_ag, het_ag, pred_ag), names_to = "modifier", values_to = "aggression")%>%
select(ind_id, obs_id, aggression, modifier)%>%
drop_na(aggression)

responses <- responses%>%
  left_join(aggression)%>%
  mutate(aggression = ifelse(is.na(aggression), 0, aggression))

# clean boris data

boris <- boris%>%
select(Observation.id, Coding.duration, Behavior,Modifier..1, Behavior.type, Duration..s., Start..s., Stop..s.)%>%
# make all colnames lower case
rename_all(tolower)%>%
# remove '.'
rename_all(~ gsub("\\.", "", .))%>%
group_by(observationid, behavior)%>%
mutate(count = n())%>%
distinct(observationid, codingduration, behavior, behaviortype, modifier1, starts, stops, durations, count)

boris_state <- boris%>%
  filter(behaviortype == "STATE")%>%
  select(observationid, codingduration,behavior, durations)%>%
  group_by(observationid, codingduration, behavior)%>%
  summarise(durations = sum(durations))%>%
  ungroup()%>%
  pivot_wider(names_from = behavior, values_from = durations)%>%
  replace(is.na(.), 0)%>%
  rename(ind_id = observationid,
    foraging_time = Foraging,
    vigilance_time = Vigilant,
    move_time = Moving)%>%
    # change types
    mutate(foraging_time = as.numeric(foraging_time),
           vigilance_time = as.numeric(vigilance_time),
           move_time = as.numeric(move_time))

boris_point = boris%>%
    filter(behaviortype == "POINT")%>%
    group_by(observationid, codingduration, behavior, modifier1)%>%
    summarise(count = last(count))%>%
    ungroup()%>%
    pivot_wider(names_from = behavior, values_from = count)%>%
    replace(is.na(.), 0)%>%
    rename(ind_id = observationid,
           bites = Bite,
           escape = Escape,
           aggression = Agression,
           modifier = modifier1)%>%
    # change types
    mutate(bites = as.numeric(bites),
           escape = as.numeric(escape),
           aggression = as.numeric(aggression))%>%
    group_by(ind_id)%>%
    summarise(bites = sum(bites),
              escape = sum(escape),
              aggression = sum(aggression),
              modifier = last(modifier),
              codingduration = max(codingduration))%>%
    ungroup()

responses_boris <- boris_state%>%
  left_join(boris_point)

# join responses

responses <- responses%>%
full_join(responses_boris)%>%
# fill all numeric columns with 0
mutate_if(is.numeric, ~replace_na(., 0))

# fix foraging time errors

errors_f_t <- responses%>%
filter(foraging_time == 0 | is.na(foraging_time) | foraging_time < 1, bites > 0)

corrected_f_t <- boris%>%
filter(observationid %in% errors_f_t$ind_id)%>%
filter(behavior == "Bite")%>%
group_by(observationid)%>%
#filter first and last bite
filter(row_number() == 1 | row_number() == n())%>%
# calculate total foraging time
reframe(foraging_time = starts - lag(starts))%>%
ungroup()%>%
filter(!is.na(foraging_time))

## add corrected foraging time

responses <- responses%>%
left_join(corrected_f_t, by = c("ind_id" = "observationid"))%>%
mutate(foraging_time = ifelse(foraging_time.x == 0 | foraging_time.x < 1, foraging_time.y, foraging_time.x))%>%
select(-foraging_time.y, -foraging_time.x)%>%
# single bites
mutate(foraging_time = ifelse(bites == 1 & is.na(foraging_time), 1, foraging_time))%>%
mutate_if(~is.numeric(.), ~replace_na(., 0))

responses <- responses%>%
select(-c(obs_id, modifier))%>%
  group_by(ind_id)%>%
  summarise(across(everything(), sum, na.rm = TRUE))%>%
  ungroup()

# remove temporary objects

rm(benthic, boris_point, boris, boris_state, category_names, corrected_f_t, damsels, errors_f_t, individuals,
   mean_rugosity, predator, rugosity, sampling, sites, aggression, behaviors, groups, plots, responses_boris, samples)

# final data

head(responses)
head(predictors)
head(fish)
head(traits)

