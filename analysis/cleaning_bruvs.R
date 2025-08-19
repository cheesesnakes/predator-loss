pacman::p_load(here, dplyr, tidyr, stringr)

# Load data

metadata <- read.csv(here("data/chapter-3_bruvs_metadata.csv"))
fish <- read.csv(here("data/chapter-3_bruvs_fish.csv"))
traits <- read.csv(here("data/fish-traits.csv"))
individuals <- read.csv(here("data/chapter-3_bruvs_size.csv"))

# guild matrix

colnames(traits) <- str_to_lower(colnames(traits))

guild <- traits %>%
    mutate_if(is.character, str_trim) %>%
    mutate(fish_spp = paste(genus, species, sep = " ")) %>%
    rename(guild = feeding.guild) %>%
    select(fish_spp, guild) %>%
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
    group_by(fish_spp) %>%
    mutate(value = value / sum(value)) %>%
    pivot_wider(names_from = guild, values_from = value) %>%
    mutate(
        Invertivore = ifelse(Invertivore < 1, 0, Invertivore),
        Piscivore = ifelse(Piscivore > 0, 1, Piscivore)
    )

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

# filter out extra species from guild data
guild <- guild %>%
filter(fish_spp %in% fish$fish_spp)%>%
select(fish_spp, Invertivore, Piscivore)

# predictors data

predictors <- metadata %>%
    select(deployment_id, location, site, protection)

# final data

head(individuals)
head(guild)
head(predictors)
