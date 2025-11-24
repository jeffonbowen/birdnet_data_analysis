
library(DBI)
library(RSQLite)
library(dplyr)

# 1. Connect
con <- dbConnect(SQLite(), "data/birds.db")

# 2. Explore
dbListTables(con)
dbListFields(con, "detections")
head(dbReadTable(con, "detections"))

# 3. Create a lazy table reference
det <- tbl(con, "detections")

# 5. Daily totals
daily <- det %>%
  group_by(Date) %>%
  summarise(n = n()) %>%
  arrange(Date) %>%
  collect()

# 6. Species totals
species_totals <- det %>%
  group_by(Com_Name) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  collect()

# 7. Save derived products for dashboard
write.csv(daily, "data/daily_totals.csv", row.names = FALSE)
write.csv(species_totals, "data/species_totals.csv", row.names = FALSE)

# 8. Disconnect
dbDisconnect(con)
