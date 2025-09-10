pacman::p_load(here, dplyr, tidyr, stringr)

# Load data

metadata <- read.csv(here("data/chapter-3_bruvs_metadata.csv"))
fish <- read.csv(here("data/chapter-3_bruvs_fish.csv"))
traits <- read.csv(here("data/fish-traits.csv"))
individuals <- read.csv(here("data/chapter-3_bruvs_size.csv"))

# guild matrix

colnames(traits) <- str_to_lower(colnames(traits))

traits %>%
    mutate_if(is.character, str_trim) %>%
    mutate(fish_spp = paste(genus, species, sep = " ")) -> traits

guild <- traits %>%
    rename(guild = feeding.guild) %>%
    select(family, fish_spp, guild) %>%
    mutate(
        FG1 = str_split(guild, ", ", simplify = T)[, 1],
        FG2 = str_split(guild, ", ", simplify = T)[, 2]
    ) %>%
    select(-c(guild)) %>%
    pivot_longer(cols = c(FG1, FG2), names_to = "type", values_to = "guild") %>%
    select(-c(type)) %>%
    mutate(guild = str_to_sentence(str_trim(guild))) %>%
    distinct() %>%
    filter(guild != "", guild != "#N/a", guild != "Plankivore", guild != "Planktivore") %>%
    mutate(value = 1) %>%
    complete(fish_spp, guild, fill = list(value = 0)) %>%
    group_by(family, fish_spp) %>%
    mutate(value = value / sum(value)) %>%
    mutate(value = ifelse(guild == "Piscivore" & value > 0, 1,
        ifelse(guild == "Invertivore" & value < 1, 0, value)
    )) %>%
    filter(value > 0)

# species with no guild information

no_guild <- fish %>%
    filter(!fish_spp %in% guild$fish_spp) %>%
    select(fish_spp) %>%
    distinct()

# remove species with no guild information from fish data
fish <- fish %>%
    filter(fish_spp %in% guild$fish_spp)

individuals <- individuals %>%
    filter(fish_spp %in% guild$fish_spp)

# Make size classes 10cm instead of 5cm, truncated at 80cm

individuals <- individuals %>%
    mutate(size = case_when(
        size %in% c("0-5", "5-10") ~ "0-10",
        size %in% c("10-15", "15-20") ~ "10-20",
        size %in% c("20-25", "25-30") ~ "20-30",
        size %in% c("30-35", "35-40") ~ "30-40",
        size %in% c("40-45", "45-50") ~ "40-50",
        size %in% c("50-55", "55-60") ~ "50-60",
        size %in% c("60-65", "65-70") ~ "60-70",
        size %in% c("70-75", "75-80") ~ "70-80",
        TRUE ~ ">80"
    ))

# filter out extra species from guild data
guild <- guild %>%
    filter(fish_spp %in% fish$fish_spp) %>%
    select(fish_spp, guild, family)

# predictors data

predictors <- metadata %>%
    select(deployment_id, location, site, protection)

# Aggregate counts by guild

individuals %>%
    left_join(guild, by = "fish_spp") %>%
    group_by(deployment_id, guild, size) %>%
    summarise(count = sum(count)) %>%
    ungroup() -> counts

# Aggregate counts by family

families <- traits %>%
    select(family, fish_spp) %>%
    distinct()

individuals %>%
    left_join(families, by = "fish_spp") %>%
    group_by(deployment_id, family, size) %>%
    summarise(count = sum(count)) %>%
    ungroup() -> family_counts

# final data

head(family_counts)
head(counts)
head(individuals)
head(guild)
head(predictors)
