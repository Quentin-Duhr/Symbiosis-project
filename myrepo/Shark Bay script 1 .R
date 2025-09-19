######################### STEP 2.1 ################################
Library(dplyr)

heatmap_data_SB.input <- read.csv("heatmap_data_SB.csv", header = TRUE)

#########################################
#heatmap_data_SB.csv = data from STAN : 
# -  prev_obs_p: initial larval abundances by genus
# - survivorship curves
# - potential encounters (from Preference_Stan.R)
# Experimental data :  day /  species (host) /  Symbiont / Temp  /  potential encounters
##########################################
#POTENTIAL OF ENCOUNTERS 
#Pe = Mean_Scaled_Preference / 1000 * Survival_Estimate
# # Reflects mosquito biting preference & observed bird abundance (aligned with wsi = global factor of infection success # # 
# Dividing by 1000 scales to meaningful range # 
# Survival_Estimate adjusts for survivability
#This is a reasonable estimate of encounter rates but does not scale across bird species to ensure it represents a proportion (e.g., the sum of wsi over all species should equal 1).#
###############
###############
#NORMALISATION OF ENCOUTER POTENTIAL TO OBTAIN wsi 
#sum of probability = 1 
###############
Scaled_potential_encounters<- heatmap_data_SB.input %>%
  group_by(Day) %>%
  mutate(Scaled_Potential_Encounters = Potential_Encounters / sum(Potential_Encounters, na.rm = TRUE))

# delta = probability infection per encounter, constant, each symbiont can only infect once
delta <- 1 
#infectious period for larvae (infectious time : 244 days max) 
D <- 244
############################################################################
############# EXTRACTION OF GLOBAL FACTOR OF INFECTION SUCCES ###############
###########################################################################

# extraction in function of different conditions (T°, Symbionte, Day) 
# FOR 27 DEGREES 
####################################################
# EX : for 27 degrees, Symbionte D and Day 7 
omega_Si_27_D.in <- Scaled_potential_encounters %>%  filter(Temp=="27", Symbiont == "D", Day == 7)
omega_Si_27_D <- omega_Si_27_D.in$Scaled_Potential_Encounters
omega_Si_27_D <- as.numeric(omega_Si_27_D)
####################################################
#2. T° = 27, Sym = C, 7 days
omega_Si_27_C.in <- Scaled_potential_encounters %>%
  filter(Temp %in% c("27")) %>% filter(Symbiont %in% c("C"))  %>% filter(Day %in% c(7)) #dim should be 16
omega_Si_27_C <- omega_Si_27_C.in$Scaled_Potential_Encounters
omega_Si_27_C <- as.numeric(omega_Si_27_C)

#3. T° = 27, Sym = D, 14 days
omega_Si_27_D.in_14 <- Scaled_potential_encounters %>%
  filter(Temp %in% c("27")) %>% filter(Symbiont %in% c("D"))  %>% filter(Day %in% c(14)) #dim should be 16
omega_Si_27_D_14 <- omega_Si_27_D.in_14$Scaled_Potential_Encounters
omega_Si_27_D_14 <- as.numeric(omega_Si_27_D_14)

#4. T° = 27, Sym = C, 14 days
omega_Si_27_C.in_14 <- Scaled_potential_encounters %>%
  filter(Temp %in% c("27")) %>% filter(Symbiont %in% c("C"))  %>% filter(Day %in% c(14)) #dim should be 16
omega_Si_27_C_14 <- omega_Si_27_C.in_14$Scaled_Potential_Encounters
omega_Si_27_C_14 <- as.numeric(omega_Si_27_C_14)

# HOT TEMP
# FOR 31 DEGREES 

#5. T° =  31, Sym = D, 7 days
omega_Si_27_D.in <- Scaled_potential_encounters %>%
  filter(Temp %in% c("31")) %>% filter(Symbiont %in% c("D"))  %>% filter(Day %in% c(7)) #dim should be 16
omega_Si_27_D <- omega_Si_31_D.in$Scaled_Potential_Encounters
omega_Si_27_D <- as.numeric(omega_Si_31_D)

#6. T° = 31, Sym = C, 7 days
omega_Si_27_C.in <- Scaled_potential_encounters %>%
  filter(Temp %in% c("31")) %>% filter(Symbiont %in% c("C"))  %>% filter(Day %in% c(7)) #dim should be 16
omega_Si_27_C <- omega_Si_27_C.in$Scaled_Potential_Encounters
omega_Si_27_C <- as.numeric(omega_Si_31_C)
########## 14 DAYS ############
#7. T° = 27, Sym =D, 14 days
omega_Si_27_D.in_14 <- Scaled_potential_encounters %>%
  filter(Temp %in% c("31")) %>% filter(Symbiont %in% c("D"))  %>% filter(Day %in% c(14)) #dim should be 16
omega_Si_27_D_14 <- omega_Si_31_D.in_14$Scaled_Potential_Encounters
omega_Si_27_D_14 <- as.numeric(omega_Si_31_D_14)

# HOT TEMP
# FOR 31 DEGREES 

#8. T° = 31, Sym = C, 14 days
omega_Si_31_C.in_14 <- Scaled_potential_encounters %>%
  filter(Temp %in% c("31")) %>% filter(Symbiont %in% c("C"))  %>% filter(Day %in% c(14)) #dim should be 16
omega_Si_31_C_14 <- omega_Si_31_C.in_14$Scaled_Potential_Encounters
omega_Si_31_C_14 <- as.numeric(omega_Si_31_C_14)
################################################################################
# HERE WE HAVE THE OMEGA OF EACH SCENARIO 
################################################################################


######################### STEP 2.2 ################################
########################################
# FOLLOW THE EVOLUTION OF PMB IN THE TIME 
########################################

# PMBd = "how many of symbiont type d are available for a larva to pick up".
# EX : INITIAL PMB calculated for GBR. 

PMB_D1a_control.dose.GBR<-1.14E+15

## Logistic functions describing the evolution of PMB over time
# Here, same formula for 27°C and 31°C, for symbionts D and C
PMB_27_D <- function(d, initial_value) {
  initial_value / (1 + exp(-0.1 * (d - 10)))
}
# => These functions represent PMBd in the equation

PMB_27_C <- function(d, initial_value) {
  initial_value / (1 + exp(-0.1 * (d - 10)))
}

PMB_31_D <- function(d, initial_value) {
  initial_value / (1 + exp(-0.1 * (d - 10)))
}

PMB_31_C <- function(d, initial_value) {
  initial_value / (1 + exp(-0.1 * (d - 10)))
}
################################################################################
# => this code defines how PMBd changes with time (d = days) instead of staying constant
################################################################################

######################### STEP 2.3 ################################
################################################################################
# SYMBIONT SURVIVAL RATES / USING CULTURED CELLS RATES AS PROXY FOR FREE-LIVING
################################################################################
######################### DATA MANAGEMENT ############################

# DATA FROM camp et al. = D1a & C1-124 Control = 26°C and Treatment = 32.4°C
#DATA FROM Chakravarti et al; Control = 27°C and Treatment = 30°C
##################
#CAMP# 
#Extracted data
#df_Camp_2022 <- data.frame(
#    Isolate = c("Dla", "Dla", "Dla", "Dla", "B1", "B1", "B1", "B1", "C1-124", "C1-124", "C1-124", "C1-124"),
#    Time_Point = c("TO", "TE", "TO", "TE", "TO", "TE", "TO", "TE", "TO", "TE", "TO", "TE"),
#    Treatment = c("Control", "Control", "Treatment", "Treatment", "Control", "Control", "Treatment", "Treatment", "Control", "Control", "Treatment", "Treatment"),
#    Cell_mL_1_Mean = c(164134, 126282, 173596, 78151, 156750, 155405, 131500, 72441, 138675, 121877, 122000, 49925),
#    Cell_mL_1_SE = c(8436, 12680, 17762, 7263, 8413, 19992, 19506, 4191, 8405, 18468, 19877, 3544),
#    F_F_Mean = c(0.55, 0.43, 0.53, 0.36, 0.53, 0.51, 0.53, 0.10, 0.59, 0.50, 0.58, 0.33),
#    F_F_SE = c(0.01, 0.01, 0.00, 0.00, 0.00, 0.01, 0.01, 0.03, 0.02, 0.00, 0.02, 0.02))

#CHAK#
## Given values from Chakravarti, need to covert growth rates to cell numbers
#initial_cells <- 200000
#growth_rate <- 0.087
#days <- 17

#initial_cells <- 200000
#growth_rate <- 0.124
#days <- 17

#initial_cells <- 200000
#growth_rate <- 0.1 
#days <- 17

#initial_cells <- 200000
#growth_rate <- 0.13
#days <- 17
# Calculate the final population
#final_population <- initial_cells * 2^(growth_rate * days)
#final_population
######################################################################
# Define inputs for each combination
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
# Function to calculate decay rate constant k,  
#############################
calculate_survival_curve <- function(T0, TE, Days) {
  # constant of decay 
  k <- -log(TE / T0) / Days
  #theorique survival function 
  survivorship_curve <- function(d) T0 * exp(-k * d)
  #equation as a text
  equation <- paste0("T(d) = ", round(T0, 2), " * exp(-", round(k, 5), " * d)")
  time_points <- seq(0, Days, by = 1)
  survival_values <- survivorship_curve(time_points)
  #Return a list (curve & equation)
  list(
    curve = data.frame(Day = time_points, Survival = survival_values),
    equation = equation
  )
}
###############################
# Generate survival curves and equations for all combinations thanks to function above
results <- lapply(inputs, function(input) {
  res <- calculate_survival_curve(input$T0, input$TE, input$Days)
  res$curve$Isolate <- input$Isolate
  res$curve$Treatment <- input$Treatment
  list(curve = res$curve, equation = res$equation)
})
###############################
# Combine curves into a single data frame
result_df <- do.call(rbind, lapply(results, function(res) res$curve))

# Create an object to store survival equation
equations <- sapply(results, function(res) res$equation)

# Display results for the survival equations 
print(result_df)    # Survival data
print(equations)    # Equations

####################################
#create a figure for the equations
# Create a data frame for survival curves
time_points <- seq(0, 244, by = 1)  # From 0 to 244 days
######
# Calculate manually the survival curve from the equations
survival_values <- list(
  Dla_Control = 164134 * exp(-0.02622 * time_points),
  Dla_Treatment = 173596 * exp(-0.07981 * time_points),
  C1_124_Control = 138675 * exp(-0.01291 * time_points),
  C1_124_Treatment = 122000 * exp(-0.08935 * time_points),
  F1_Control = 200000 * exp(-0.09011 * time_points),
  F1_Treatment = 200000 * exp(-0.06931 * time_points),
  G3_Control = 200000 * exp(-0.08595 * time_points),
  G3_Treatment = 200000 * exp(-0.0603 * time_points)
)
######################################
# Convert the list into a data frame
survival_df <- data.frame(
  Time = rep(time_points, length(survival_values)),
  Survival = unlist(survival_values),
  Isolate = rep(c("Dla", "Dla", "C1-124", "C1-124", "F1", "F1", "G3", "G3"), each = length(time_points)),
  Treatment = rep(c("Control", "Treatment"), times = 4, each = length(time_points))
)

######################## Plot the survival curves############################
ggplot(survival_df, aes(x = Time, y = Survival, color = interaction(Isolate, Treatment))) +
  geom_line(size = 1) +
  labs(title = "Survivorship Curves",
       x = "Time (Days)",
       y = "Population Size",
       color = "Model (Isolate, Treatment)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "purple", "green", "lightblue", "red", "red", "pink", "yellow"))

#######################################################################
#####################DEFINE SURVIVAL FUNCTION##########################
#######################################################################
#2. free-living survival curves, specific per symbiont taxa, derived from experiments with free-living/cultured cells
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

################################################################################
#################################################################
######################### STEP 2.5 ################################
################################################################################
#at this point, it is only set up for 27C and D at day 7
# Total de symbionts D (mosquitos), à 27°C, calculé pour chaque jour (1:D).
# Ce total combine :
#  - la fonction de ponte/production (PMB_27_D),
#  - la dose initiale de symbionts dans l’environnement (PMB_D1a_control.dose.GBR),
#  - la survie des symbionts au cours du temps (SM_27_D1a)total_mosquitoes <- sapply(1:D, 
function(d) {
  PMB_27_D(d, PMB_D1a_control.dose.GBR) * SM_27_D1a(d)
}


###########################INFECTIOUS POTENTIAL PER SPECIES##################
## EQUATION 1 – PART 1 
mu_i_27_D.test4_totM <- sapply(1:length(omega_Si_27_D), function(i) {  
  # Infectious potential of the species (i)
  potential_infected <- sapply(1:D, function(d) {
    # Random effect (uniform) = stochastic variability 
    random_factor <- runif(1, 0.5, 1.5)
    # Infectious potential
    omega_Si_27_D[i] * pref_abund_subset$Potential_Encounters[i] * random_factor
  })
  # # replace NA by 0
  potential_infected[is.na(potential_infected)] <- 0
  # # Cap the total potential infected for species i
  scaling_factor <- min(1, pref_abund_subset$Potential_Encounters[i] / sum(potential_infected, na.rm = TRUE))
  # Scale the daily potentials to ensure they respect the cap
  potential_infected <- potential_infected * scaling_factor
  # Ensure the column sum doesn't exceed the species' total potential encounters
  while (sum(potential_infected, na.rm = TRUE) > pref_abund_subset$Potential_Encounters[i]) {
    potential_infected <- potential_infected * (pref_abund_subset$Potential_Encounters[i] / sum(potential_infected, na.rm = TRUE))
  }
  return(potential_infected)
})
###########################DAY TRANSMISSION EFFECTIVE######################
## EQUATION 1 – PART 2
daily_transmission_potential <- apply(mu_i_27_D.test4_totM, 2, function(day_totals) {
  #scaling in function of the number of mosquitoes
  scaling_factor <- min(1, total_mosquitoes / sum(day_totals, na.rm = TRUE))
  #factor application
  day_totals * scaling_factor
})
#PLOT TO SEE IF THAT MAKE SENSE 
#look at daily infection ? variable by day? variable by species?
#Convert matric in DF #
mu_i_27_D.test4_totM.df<-as.data.frame(daily_transmission_potential)
#Add name of species
species_names <- pref_abund_subset$Species
colnames(mu_i_27_D.test4_totM.df) <- species_names
#ADD column day
mu_i_27_D.test4_totM.df <- mutate(mu_i_27_D.test4_totM.df, Day = row_number())

# Convert data to long format for ggplot
mu_i_27_D.test4_totM.df.reshape <- pivot_longer(mu_i_27_D.test4_totM.df, cols = -c("Day"),  # Exclude "Day" from reshaping
                                                names_to = "Species", values_to = "Daily_susceptible_and_infected_birds") 

######## VISUALISATION (quotidian valour per species) 
ggplot(mu_i_27_D.test4_totM.df.reshape, aes(x = Day, y = Daily_susceptible_and_infected_birds, color = Species)) + 
  geom_point() + geom_line()+
  facet_wrap(~Species, scales = "free_y") + 
  labs(y = "Count", color = "Species", title ="Daily number of infected larvae via transmission") + 
  theme_minimal()

################################################################################ 
#TOTAL INFECTION IN THE TIME PER SPECIES
################################################################################ 

## Cumulative values of infected birds do not exceed total available birds
mu_i_27_D.test4_cumulative <- apply(mu_i_27_D.test4_totM, 2, cumsum)
#check that these values make sense by plotting
# Create a data frame combining the mosquito data and the bird infection data
# Convert matrix to data frame
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
#DATA FUSION FOR CAMPARAISON 
################################################################################ 
# prepare final df for comparison
###########
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
