pacman::p_load(here, dplyr, tidyr, stringr)
#
# bruvs <- read.csv(here("data/chapter-3_bruvs.csv"))
#
# head(bruvs)
# colnames(bruvs)
#
# metadata <- bruvs %>%
#     select(Date, deployments, region, Site, status) %>%
#     mutate_if(is.character, str_trim) %>%
#     distinct() %>%
#     rename(deployment_id = deployments)
#
# colnames(metadata) <- tolower(colnames(metadata))
#
# head(metadata)
#
# write.csv(metadata, here("data/chapter-3_bruvs_metadata.csv"), row.names = FALSE)
#
# fish <- bruvs %>%
#     select(deployments, Fish.spp, Media.Code, T1st, MaxN, MaxN.Time, MxMC, total_fish) %>%
#     mutate_if(is.character, str_trim) %>%
#     rename(
#         deployment_id = deployments,
#         fish_spp = Fish.spp,
#         media_code = Media.Code,
#         time_first = T1st,
#         max_n = MaxN,
#         max_n_time = MaxN.Time,
#     )
#
# colnames(fish) <- tolower(colnames(fish))
#
# head(fish)
#
# write.csv(fish, here("data/chapter-3_bruvs_fish.csv"), row.names = FALSE)
#
# traits <- bruvs %>%
#     select(genus, Fish.spp, guild) %>%
#     mutate_if(is.character, str_trim) %>%
#     distinct() %>%
#     rename(
#         fish_spp = Fish.spp,
#         genus = genus,
#         guild = guild
#     )
#
# colnames(traits) <- tolower(colnames(traits))
# head(traits)
#
# write.csv(traits, here("data/chapter-3_bruvs_traits.csv"), row.names = FALSE)
#
# size <- bruvs %>%
#     select(-c(Date, region, Site, status, genus, guild, Media.Code, T1st, MaxN, MaxN.Time, total_fish, MxMC)) %>%
#     select(-starts_with("X")) %>%
#     rename(
#         deployment_id = deployments,
#         fish_spp = Fish.spp
#     ) %>%
#     pivot_longer(
#         cols = starts_with("Fish."),
#         names_to = "Fish",
#         values_to = "size"
#     ) %>%
#     select(-Fish) %>%
#     filter(!is.na(size) & size != "" & size != " ") %>%
#     group_by(deployment_id, fish_spp, size) %>%
#     summarise(
#         count = n(),
#         .groups = "drop"
#     ) %>%
#     mutate_if(is.character, str_trim)
#
# colnames(size) <- tolower(colnames(size))
# head(size)
#
# write.csv(size, here("data/chapter-3_bruvs_size.csv"), row.names = FALSE)
#

# Load data

metadata <- read.csv(here("data/chapter-3_bruvs_metadata.csv"))
fish <- read.csv(here("data/chapter-3_bruvs_fish.csv"))
traits <- read.csv(here("data/chapter-3_bruvs_traits.csv"))
individuals <- read.csv(here("data/chapter-3_bruvs_size.csv"))

# guild matrix

guild <- traits %>%
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

# predictors data

predictors <- metadata %>%
    select(deployment_id, region, site, status) %>%
    rename(
        location = region,
        protection = status
    )

# final data

head(individuals)
head(guild)
head(predictors)
