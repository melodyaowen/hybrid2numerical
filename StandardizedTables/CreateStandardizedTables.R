source("./RequiredPackages.R")

create_standardized_table <- function(file_in = NA,
                                      file_out = NA){

  myDat <- read.csv(file_in) %>%
    dplyr::select(-group_id, -eff1eff2, -rho01rho02) %>%
    rename(group_id = "X") %>%
    mutate(across(contains("method"),
                  ~ case_when(str_detect(.x, "n = ") ~ as.numeric(str_extract(.x, "(?<=n = )\\d+")),
                              .x == "0%" ~ 0,
                              TRUE ~ NA_real_)))

  tableDat <- myDat %>%
    #filter(eff2minus1 >= 0) %>%
    mutate(
      eff2minus1_group = case_when(
        eff2minus1 < 0 ~ "-",
        eff2minus1 == 0 ~ "0",
        eff2minus1 > 0    & eff2minus1 <= 0.19 ~ "0.05 to 0.19",
        eff2minus1 > 0.19 & eff2minus1 <= 0.29 ~ "0.20 to 0.29",
        eff2minus1 > 0.29 & eff2minus1 <= 0.39 ~ "0.30 to 0.39",
        eff2minus1 > 0.39 & eff2minus1 <= 0.49 ~ "0.40 to 0.49",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      rho02minus01_group = case_when(
        rho02minus01 == 0 ~ "rho01 = rho02",
        rho02minus01 > 0  ~ "rho01 < rho02",
        rho02minus01 < 0  ~ "rho01 > rho02",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(eff2minus1_group, rho02minus01) %>%
    dplyr::select(-group_id, -eff2minus1, -rho02minus01) %>%
    relocate(eff2minus1_group, rho02minus01_group) %>%
    group_by(eff2minus1_group, rho02minus01_group) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(group_id = row_number()) %>%
    relocate(group_id) %>%
    group_by(group_id) %>%
    mutate(across(contains("method"),
                  ~ if_else(.x == 0, "0%",  # the special case: just "0%"
                            paste0(round(.x/n*100, 0), "% (n = ", .x, ")"))))

  write.csv(tableDat,
            file = paste0("./StandardizedTables/", file_out))

}

# Comparison I
create_standardized_table(file_in = "Comparison1/Results1/BestAll_1_STD.csv",
                          file_out = "Summary_1_STD.csv")

# Comparison II
create_standardized_table(file_in = "Comparison2/Results2/BestAll_2_STD.csv",
                          file_out = "Summary_2_STD.csv")

# Comparison III
create_standardized_table(file_in = "Comparison3/Results3/BestAll_3_STD.csv",
                          file_out = "Summary_3_STD.csv")

# Comparison IV
create_standardized_table(file_in = "Comparison4/Results4/BestAll_4_STD.csv",
                          file_out = "Summary_4_STD.csv")
