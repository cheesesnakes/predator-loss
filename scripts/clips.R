setwd(".")

# calculating location of sample in its respective clips

samples <- read.csv("../data/chapter-3_samples.csv")

head(samples)

# set length of videos at 11:47 in seconds

l <- 11*60 + 47

l

# make frame of start and end of 8 clips in steps of l

clips <- data.frame(clip_id = 1:8, start = seq(0, l*7, by = l), end = seq(l, l*8, by = l))


# use sample_s to locate sample in clips and add clip_id
library(dplyr)
library(tidyr)

samples <- samples %>%
  mutate(clip_id = case_when(
    sample_s <= l ~ 1,
    sample_s <= 2*l ~ 2,
    sample_s <= 3*l ~ 3,
    sample_s <= 4*l ~ 4,
    sample_s <= 5*l ~ 5,
    sample_s <= 6*l ~ 6,
    sample_s <= 7*l ~ 7,
    sample_s <= 8*l ~ 8
  ))

# add second in clip

samples <- samples %>%
  left_join(clips, by = "clip_id") %>%
  mutate(second_in_clip = sample_s - start)

# convert to minutes and seconds

samples <- samples %>%
  mutate(minute = floor(second_in_clip/60),
         second = second_in_clip %% 60)

# join minute and second to clip_time

samples <- samples %>%
  mutate(clip_time = paste0(minute, ":", second))

samples <- samples%>%
    dplyr::select(deployment_id, plot_id, sample_id, sample_s, sample_m, clip_id, clip_time, completed, remarks)

head(samples)

write.csv(samples, "chapter-3_samples.csv", row.names = FALSE)
