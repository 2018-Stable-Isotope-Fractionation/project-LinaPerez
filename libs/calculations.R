# saturation vapor pressure in [Pa] from http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html
# this is the one we use in the model
calc_es <- function(Temp.K) {
  exp(77.3450 + 0.0057 * Temp.K - 7235 / Temp.K) / Temp.K^8.2  
}

# saturation vapor pressure in [Pa], 
# slightly different equation from http://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/qsat.html
# not used in the model
calc_es2 <- function(Temp.K) {
  Temp.C <- calc_Temp.C(Temp.K)
  (6.1121*(1.0007) * exp((17.502*Temp.C) / (240.97+Temp.C))) * 100 # mb to Pa
}

# derivative of saturation vapor pressure
calc_des_dT <- function(Temp.K) {
  calc_es(Temp.K) * (0.0057 + 7235 / Temp.K^2) - 8.2 * calc_es(Temp.K) / Temp.K
}

# temperature conversion convenience functions
calc_Temp.K <- function(Temp.C) 273.15 + Temp.C
calc_Temp.C <- function(Temp.K) Temp.K - 273.15

# calculate the latent heat of water condensation (polynomial valid from -25 to 40) in [J/kg]
calc_L <- function(Temp.K) {
  Temp.C <- calc_Temp.C (Temp.K)
  (2500.8 - 2.36 * Temp.C + 0.0016 * Temp.C^2 - 0.00006*Temp.C^3 )*1000 # in J/kg
}

# calculate saturation mass mixing ratio of water i.e. the saturation fraction of vapor in the total gas
# also called saturation vapor fraction (q = vapor fraction)
# es is saturation wapour pressure in Pa, p is pressure in Pa
# return dimensionless (i.e. fraction) ratio
calc_qs <- function(es, p) 0.622 * es/p
calc_q_from_RH <- function(RH, es, p) RH * calc_qs(es, p)

# fractionation factors
calc_O_alpha_l_v <- function(Temp.K) exp( (-7.685 + 6.7123 * 10^3/Temp.K - 1.6664 * 10^6/Temp.K^2 + 0.35041 * 10^9 / Temp.K^3) / 1000 )
calc_O_alpha_i_v <- function(Temp.K) exp( (11839/Temp.K - 28.224)/1000 )
calc_O_alpha <- function(Temp.K) {
  fraction_water <- (Temp.K - 253.15)/(273.15 - 253.15)
  fraction_water <- ifelse(fraction_water < 0, 0, ifelse(fraction_water > 1, 1, fraction_water))
  fraction_water * calc_O_alpha_l_v (Temp.K) + (1-fraction_water) * calc_O_alpha_i_v (Temp.K)
}

# for deuterium
calc_H_alpha_l_v <- function(Temp.K) exp( (1158.8 *Temp.K^3/ 10^9 - 1620.1 *Temp.K^2/ 10^6 + 794.84 * Temp.K/ 10^3 - 161.09 + 2.9992*10^9/Temp.K^3) / 1000 )
calc_H_alpha_i_v <- function(Temp.K) exp( (16288/Temp.K^2 *10^3 - 93.4)/1000 )
calc_H_alpha <- function(Temp.K) {
  fraction_water <- (Temp.K - 253.15)/(273.15 - 253.15)
  fraction_water <- ifelse(fraction_water < 0, 0, ifelse(fraction_water > 1, 1, fraction_water))
  fraction_water * calc_H_alpha_l_v (Temp.K) + (1-fraction_water) * calc_H_alpha_i_v (Temp.K)
}

