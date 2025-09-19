      ######################### STEP 2.1 ################################
library(dplyr)
library(ggplot2)
heatmap_data_SB.input <- read.csv("heatmap_data_SB.csv", header = TRUE)

################ DATA SUMMARY #############

#heatmap_data_SB.csv = data from STAN : 
# - prev_obs_p: initial larval abundances by genus
# - survivorship curves
# - potential encounters (from Preference_Stan.R)
# Experimental data :  day /  species (host) /  Symbiont / Temp  /  potential encounters

#data#########################################

#POTENTIAL OF ENCOUNTERS 
#Pe = Mean_Scaled_Preference / 1000 * Survival_Estimate
# # Reflects mosquito biting preference & observed bird abundance (aligned with wsi = global factor of infection success # # 
# Dividing by 1000 scales to meaningful range # 
# Survival_Estimate adjusts for survivability
#This is a reasonable estimate of encounter rates but does not scale across bird species to ensure it represents a proportion (e.g., the sum of wsi over all species should equal 1).#


################ NORMALISATION  #############

# Normalize encounter potentials to obtain ωSi (sum = 1 per group)
Scaled_potential_encounters <- heatmap_data_SB.input %>%
  group_by(Day) %>%
  mutate(Scaled_Potential_Encounters = Potential_Encounters / sum(Potential_Encounters, na.rm = TRUE))

# delta = probability infection per encounter, constant, each symbiont can only infect once
delta <- 1 
#infectious period for larvae (infectious time : 244 days max) 
D <- 244

 #
#+# delta, D and Scaled_potential_Encounters created (Probability between 0 & 1)
 #

###########################################################################
############# EXTRACTION OF GLOBAL FACTOR OF INFECTION SUCCES ############# 
###############################(**wsi**)###################################

          #################################################
# FUNCTION to extract ωSi for any (Temp, Symbiont, Day) combination
          #################################################
get_omega <- function(data, temp, symbiont, day) {
  val <- data %>%
    filter(Temp == temp, Symbiont == symbiont, Day == day) %>%
    pull(Scaled_Potential_Encounters)
  return(as.numeric(val))
}

                          ##############
# Generate ωSi values for all combinations of Temp × Symbiont × Day
                          ##############

temps <- c("27", "31")
symbionts <- c("C", "D")
days <- c(7, 14)

for (temp in temps) {
  for (sym in symbionts) {
    for (day in days) {
      obj_name <- paste0("omega_Si_", temp, "_", sym, "_", day)
      assign(obj_name, get_omega(Scaled_potential_encounters, temp, sym, day))
    }
  }
}

################################################################################
# CALL THE OMEGA OF EACH SCENARIO 
################################################################################
print(omega_Si_27_D_7)
print(omega_Si_31_C_14)


      ######################### STEP 2.2 ################################
            ########################################
            # FOLLOW THE EVOLUTION OF PMB IN THE TIME 
            #########Logistic growth model##########

# PMBa = "how many of symbiont type a are available for a larva to pick up".

# -------------------------------
# Initial value (example GBR control dose)
# -------------------------------
PMB_D1a_control.dose.GBR <- 1.14E+15

# ================================
# General logistic PMB function
# ================================

PMB_logistic <- function(d, initial_value, k = 0.1, t0 = 10) {
  # Logistic growth/decay function
  initial_value / (1 + exp(-k * (d - t0)))
}

# ================================
# Wrapper functions for each condition
# ================================

# 27°C
PMB_27_D <- function(d, initial_value = PMB_D1a_control.dose.GBR) {
  PMB_logistic(d, initial_value, k = 0.1, t0 = 10)
}

PMB_27_C <- function(d, initial_value = PMB_D1a_control.dose.GBR) {
  PMB_logistic(d, initial_value, k = 0.1, t0 = 10)
}

# 31°C
PMB_31_D <- function(d, initial_value = PMB_D1a_control.dose.GBR) {
  PMB_logistic(d, initial_value, k = 0.1, t0 = 10)
}

PMB_31_C <- function(d, initial_value = PMB_D1a_control.dose.GBR) {
  PMB_logistic(d, initial_value, k = 0.1, t0 = 10)
}

################################################################################
# => this code defines how PMBd changes with time (d = days) 
################################################################################


      ######################### STEP 2.3 ################################

################################################################################
# SYMBIONT SURVIVAL RATES / USING CULTURED CELLS RATES AS PROXY FOR FREE-LIVING
################################################################################

    ######################### DATA MANAGEMENT ############################


# ================================
# CAMP et al. (2022) 
#D1a & C1-124 Control = 26°C and Treatment = 32.4°C
# ================================
df_Camp_2022 <- data.frame(
  Isolate = c("D1a", "D1a", "D1a", "D1a", 
              "B1", "B1", "B1", "B1", 
              "C1-124", "C1-124", "C1-124", "C1-124"),
  Time_Point = rep(c("T0", "TE"), 6),
  Treatment = rep(c("Control", "Control", "Treatment", "Treatment"), 3),
  Cell_mL_1_Mean = c(164134, 126282, 173596, 78151, 
                     156750, 155405, 131500, 72441, 
                     138675, 121877, 122000, 49925),
  Cell_mL_1_SE   = c(8436, 12680, 17762, 7263, 
                     8413, 19992, 19506, 4191, 
                     8405, 18468, 19877, 3544),
  F_F_Mean = c(0.55, 0.43, 0.53, 0.36, 
               0.53, 0.51, 0.53, 0.10, 
               0.59, 0.50, 0.58, 0.33),
  F_F_SE   = c(0.01, 0.01, 0.00, 0.00, 
               0.00, 0.01, 0.01, 0.03, 
               0.02, 0.00, 0.02, 0.02)
)


# ================================
# CHAKRAVARTI et al. 20..
#Control = 27°C and Treatment = 30°C
# ================================

# Common settings
initial_cells <- 200000
days <- 17

# ================================
# Generic growth function
# ================================
symbiont_growth <- function(initial_cells, growth_rate, days) {
  initial_cells * 2^(growth_rate * days)
}
# Define separate objects for each case 
chak_control_low   <- symbiont_growth(initial_cells, growth_rate = 0.087, days = days)
chak_control_high  <- symbiont_growth(initial_cells, growth_rate = 0.124, days = days)
chak_treatment_low <- symbiont_growth(initial_cells, growth_rate = 0.100, days = days)
chak_treatment_high<- symbiont_growth(initial_cells, growth_rate = 0.130, days = days)

      ##########################################################

# ================================
# Define raw data (T0, TE extracted or derived)
# ================================

inputs <- list(
  list(Isolate = "Dla", Treatment = "Control", T0 = 164134, TE = 126282, Days = 10),
  list(Isolate = "Dla", Treatment = "Treatment", T0 = 173596, TE = 78151, Days = 10),
  list(Isolate = "C1-124", Treatment = "Control", T0 = 138675, TE = 121877, Days = 10),
  list(Isolate = "C1-124", Treatment = "Treatment", T0 = 122000, TE = 49925, Days = 10),
  list(Isolate = "F1", Treatment = "Control", T0 = 200000, TE = 925350.5, Days = 17),
  list(Isolate = "F1", Treatment = "Treatment", T0 = 200000, TE = 649801.9, Days = 17),
  list(Isolate = "G3", Treatment = "Control", T0 = 200000, TE = 862186.5, Days = 17),
  list(Isolate = "G3", Treatment = "Treatment", T0 = 200000, TE = 557510.9, Days = 17)
)
#########################################################################

                  #############################
          #Function to calculate decay rate constant k,  
                  #############################

calculate_survival_curve <- function(T0, TE, Days) {
  #---------------------
  # Constant of decay k
  k <- -log(TE / T0) / Days
  #---------------------
  # Theoretical survival function
  survivorship_curve <- function(d) T0 * exp(-k * d)
  #---------------------
  # Equation as a text
  equation <- paste0("T(d) = ", round(T0, 2), " * exp(-", round(k, 5), " * d)")
  time_points <- seq(0, Days, by = 1)
  survival_values <- survivorship_curve(time_points)

  #------------------------------
  # Return a list (curve, equation, decay rate)
  list(
    curve = data.frame(Day = time_points, Survival = survival_values),
    equation = equation
  )
}


###############################
# Generate survival curves and equations for all combinations
# using the calculate_survival_curve function
###############################

results <- lapply(inputs, function(input) {
  #---------------------------
  # Calculate survival curve for the current input
  res <- calculate_survival_curve(input$T0, input$TE, input$Days)
  # Add metadata: Isolate and Treatment
  res$curve$Isolate <- input$Isolate
  res$curve$Treatment <- input$Treatment
  list(curve = res$curve, equation = res$equation)
} )


# Combine all curves into a single data frame
#========================
result_df <- do.call(rbind, lapply(results, function(res) res$curve))

# Extract equations for each isolate-treatment combination
#========================
equations <- sapply(results, function(res) res$equation)

# Display results for the survival equations 
print(result_df)    # Survival data
print(equations)    # Equations


####################################
# Example: manually calculated survival curves from equations
####################################

# Time sequence for manual calculation (0 to 244 days)
time_points <- seq(0, 244, by = 1)

# Survival values for each isolate-treatment combination
survival_values <- list(
  Dla_Control      = 164134 * exp(-0.02622 * time_points),
  Dla_Treatment    = 173596 * exp(-0.07981 * time_points),
  C1_124_Control   = 138675 * exp(-0.01291 * time_points),
  C1_124_Treatment = 122000 * exp(-0.08935 * time_points),
  F1_Control       = 200000 * exp(-0.09011 * time_points),
  F1_Treatment     = 200000 * exp(-0.06931 * time_points),
  G3_Control       = 200000 * exp(-0.08595 * time_points),
  G3_Treatment     = 200000 * exp(-0.06030 * time_points)
)
######################################
# Convert the list into a data frame
survival_df <- data.frame(
  Time      = rep(time_points, length(survival_values)),       # time
  Survival  = unlist(survival_values),                        # survival values
  Isolate   = rep(c("Dla", "Dla", "C1-124", "C1-124", "F1", "F1", "G3", "G3"), 
                  each = length(time_points)),                # isolate ID
  Treatment = rep(c("Control", "Treatment"), times = 4, 
                  each = length(time_points))                 # treatment type
)

######################## 
# Plot the survival curves
ggplot(survival_df, aes(x = Time, y = Survival, color = interaction(Isolate, Treatment))) +
  geom_line(size = 1) +
  labs(title = "Survivorship Curves",
       x = "Time (Days)",
       y = "Population Size",
       color = "Model (Isolate, Treatment)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "purple", "green", "lightblue", 
                                "red", "orange", "pink", "yellow"))

#######################################################################
#####################DEFINE SURVIVAL FUNCTION##########################
#######################################################################

# free-living survival curves, specific per symbiont taxa, derived from experiments with free-living/cultured cells

###############27°C#############
SM_27_D1a <- function(d) { 164134 * exp(-0.02622 * d) }
SM_C1_124_Control <- function(d) { 138675 * exp(-0.01291 * d) }
SM_F1_Control <- function(d) { 200000 * exp(-0.09011 * d) }
SM_G3_Control <- function(d) { 200000 * exp(-0.08595 * d) }

###############32°C#############
SM_Dla_Treatment <- function(d) { 173596 * exp(-0.07981 * d) }
SM_C1_124_Treatment <- function(d) { 122000 * exp(-0.08935 * d) }
SM_F1_Treatment <- function(d) { 200000 * exp(-0.06931 * d) }
SM_G3_Treatment <- function(d) { 200000 * exp(-0.0603 * d) }

######################### STEP 2.4 ################################
################################################################################
#samp_data_with_fecunditymean.SB_present_1.5_repo created in Github script: Fecundity_SB.R
#creating a DF
pref_abund_subset <- data.frame(
  Species = SCIENTIFIC.NAME$samp_data_with_fecunditymean.SB_present_1.5_repo, 
  Potential_Encounters = Total_repo_output_SB$samp_data_with_fecunditymean.SB_present_1.5_repo,
  scenario = scenario$samp_data_with_fecunditymean.SB_present_1.5_repo, 
  source = source$samp_data_with_fecunditymean.SB_present_1.5_repo)

######################### STEP 2.5 ################################


              ######################################
              # TOTAL MOSQUITOES (symbionts D at 27°C)
              ######################################
# Total number of symbionts D (mosquitos), at 27°C, for each day (1:D).
# This combines :
#  - egg laying/production function (PMB_27_D),
#  - initial dose of symbionts in environment (PMB_D1a_control.dose.GBR),
#  - survival of symbionts over time (SM_27_D1a)

total_mosquitoes <- sapply(1:D, function(d) {
  PMB_27_D(d, PMB_D1a_control.dose.GBR) * SM_27_D1a(d)
})

              ######################################
              # INFECTIOUS POTENTIAL PER SPECIES
              # EQUATION 1 – PART 1
              ######################################
mu_i_27_D.test4_totM <- sapply(1:length(omega_Si_27_D), function(i) {  
  # Infectious potential of species i
  potential_infected <- sapply(1:D, function(d) {
    # Random effect (uniform distribution) = stochastic variability by day
    random_factor <- runif(1, 0.5, 1.5)
    # Infectious potential for species i on day d
    omega_Si_27_D[i] * pref_abund_subset$Potential_Encounters[i] * random_factor
    
  })
  
  # Replace NA with 0
  potential_infected[is.na(potential_infected)] <- 0
  # Cap the total infectious potential to the species' max encounters
  scaling_factor <- min(1, pref_abund_subset$Potential_Encounters[i] / sum(potential_infected, na.rm = TRUE))
  potential_infected <- potential_infected * scaling_factor
  # Ensure the column sum never exceeds species' encounters (safety loop)
  while (sum(potential_infected, na.rm = TRUE) > pref_abund_subset$Potential_Encounters[i]) {
    potential_infected <- potential_infected * 
      (pref_abund_subset$Potential_Encounters[i] / sum(potential_infected, na.rm = TRUE))
  }
  
  return(potential_infected)
})

                ######################################
                # DAY TRANSMISSION EFFECTIVE
                # EQUATION 1 – PART 2
                ######################################
daily_transmission_potential <- apply(mu_i_27_D.test4_totM, 2, function(day_totals) {
  # Scaling by the number of mosquitoes
  scaling_factor <- min(1, total_mosquitoes / sum(day_totals, na.rm = TRUE))
  
  # Apply scaling factor
  day_totals * scaling_factor
})

                ######################################
                # DATA PREPARATION FOR PLOTTING
                ######################################

mu_i_27_D.test4_totM.df <- as.data.frame(daily_transmission_potential)

species_names <- pref_abund_subset$Species

colnames(mu_i_27_D.test4_totM.df) <- species_names

mu_i_27_D.test4_totM.df <- mutate(mu_i_27_D.test4_totM.df, Day = row_number())

# Reshape into long format for ggplot
mu_i_27_D.test4_totM.df.reshape <- pivot_longer(
  mu_i_27_D.test4_totM.df,
  cols = -c("Day"),                          
  names_to = "Species",
  values_to = "Daily_susceptible_and_infected_birds"
)

######## VISUALISATION (quotidian valour per species) #########
ggplot(mu_i_27_D.test4_totM.df.reshape, aes(x = Day, y = Daily_susceptible_and_infected_birds, color = Species)) + 
  geom_point() + geom_line()+
  facet_wrap(~Species, scales = "free_y") + 
  labs(y = "Count", color = "Species", title ="Daily number of infected larvae via transmission") + 
  theme_minimal()

################################################################################ 
################TOTAL INFECTION IN THE TIME PER SPECIES#########################
################################################################################ 

## Cumulative values of infected birds do not exceed total available birds
mu_i_27_D.test4_cumulative <- apply(mu_i_27_D.test4_totM, 2, cumsum)

mu_i_27_D.test4_cumulatived.f <- as.data.frame(mu_i_27_D.test4_cumulative)
# Get species names from omega_Si_27_D.in
species_names <- omega_Si_27_D.in$species 

# Add species names as column names to the data frame
colnames(mu_i_27_D.test4_cumulatived.f) <- species_names

# Create a new column named "Day" with row numbers
mu_i_27_D.test4_cumulatived.f <- mutate(mu_i_27_D.test4_cumulatived.f, Day = row_number())

# Convert data to long format
check_cum <- pivot_longer(mu_i_27_D.test4_cumulatived.f, cols = -c("Day"),  # Exclude "Day" from reshaping
                          names_to = "Species", values_to = "Cumulative_Value") 
# Create the ggplot
ggplot(check_cum, aes(x = Day, y = Cumulative_Value)) +
  geom_line() + facet_wrap(~Species, scales = "free_y")+scale_y_log10()+
  labs(x = "Day", y = "Cumulative Value (log 10)", title = "Cumulative Values compared to pref_abund_subset$Potential_Encounters") +
  theme_bw()+geom_hline(data = pref_abund_subset, aes(yintercept = Potential_Encounters), linetype = "dashed")

################################################################################ 
##################### DATA FUSION FOR CAMPARAISON  #############################
################################################################################ 
#=====================
# prepare final df for comparison
#=====================

data_for_plot_final <- data.frame(
  day = rep(1:D, each = length(omega_Si_27_D.in$species)),
  Species = rep(omega_Si_27_D.in$species, times = D),
  total_symbionts_D = rep(total_mosquitoes, each = length(omega_Si_27_D.in$species)),
  Total_larvae = rep(pref_abund_subset$Potential_Encounters, times = D), #Availability of potential encounters, akak, this is the number of available larvae after accounting for larval survival rates over time
  Daily_susceptible_and_infected_birds = mu_i_27_D.test4_totM.df.reshape$Daily_susceptible_and_infected_birds,
  CumTotal_suscep_and_infect_birds = check_cum$Cumulative_Value, times = D)

###################### 
# Final visualisation : dynamics comparaison
######################

ggplot(data_for_plot_final, aes(x = day)) +
  geom_line(aes(y = total_symbionts_D, color = "Total pool of available symbionts D at 27 (GBR)")) +
  geom_line(aes(y = Daily_susceptible_and_infected_birds, color = "Daily number of susceptible and infected larvae")) +
  geom_line(aes(y = Total_larvae, color = "Total pool of number of available larvae from spawning event X")) +
  geom_line(aes(y = CumTotal_suscep_and_infect_birds, color ="Cumulative total number of susceptible and infected larvae"))+
  facet_wrap(~Species, scales = "free_y")+
  scale_y_log10() +
  labs(y = "Count (log10)", color = "Legend") +
  theme_minimal()

#data_for_plot_final_SBtest.pdf  

###############################################################################
