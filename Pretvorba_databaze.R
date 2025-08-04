library(readr)
library(dplyr)
library(tidyr)

df <- read_delim("dataset.csv", delim = ";", locale = locale(encoding = "latin1"))

# Identify non-numeric columns
non_numeric_cols <- df %>% select(where(is.character)) %>% colnames()
print(non_numeric_cols)

df_long <- df %>%
  pivot_longer(
    cols = -c(Elevation, Mountain, Temperature, `Time period`, Wind),  # Keep environmental variables
    names_to = "Species",
    values_to = "Number"
  )

df_long <- df_long %>%
  filter(Number > 0)

glimpse(df_long)

summary(df_long$Number)

write_csv(df_long, "refined_dataset.csv")

######
# Load refined dataset (species occurrences)
df_occurrences <- read_csv("refined_dataset.csv")

df_traits <- read_delim("functional.csv", delim = ";", locale = locale(encoding = "latin1"))

# Clean 'Species' column in both datasets
df_occurrences <- df_occurrences %>%
  mutate(Species = str_trim(Species))

df_traits <- df_traits %>%
  mutate(Species = str_trim(Species))  # Remove leading/trailing spaces

df_final <- df_occurrences %>%
  left_join(df_traits, by = "Species")

anti_join(df_occurrences, df_traits, by = "Species")  # Show species without a match

# Save the cleaned dataset
write_csv(df_final, "final_dataset_with_traits.csv")

