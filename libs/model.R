# model base plot 
model_base_plot <- 
  ggplot() +
  aes(x = value, y = z, color= variable, linetype = scenario) +
  scale_y_continuous(breaks = seq(0, 10000, by = 1000), expand = c(0,0)) +
  geom_hline(data = function(df) {
    df %>% group_by(scenario) %>% 
      summarize(yintercept = max(z[is.na(R_c)]))
  }, aes(yintercept = yintercept, linetype = scenario), color = "black") +
  geom_line() +
  facet_wrap(~panel, scales = "free_x") + 
  theme_bw()

# run the model
run_model <- function (parameters, isotope) {
  
  # isotope system
  if (isotope == "18O") {
    ref_ratio <- get_standard(isotope) %>% get_value()
    assign("calc_alpha", calc_O_alpha, envir = .GlobalEnv)
  } else if (isotope == "2H") {
    ref_ratio <- get_standard(isotope) %>% get_value()
    assign("calc_alpha", calc_H_alpha, envir = .GlobalEnv)
  } else
    stop("don't know element: ", element)
  
  parameters %>% 
    # group by scenario
    group_by(scenario) %>% 
    # do the calculation for each parameter set
    do({
      sdf <- .
      # initial variable
      variables_initial <- data_frame(
        z = 0, # the height [m]
        G = 0, # dimensionless
        Temp.K = 273.15 + sdf$T_start[1], # [K]
        es = calc_es(Temp.K),
        q = calc_q_from_RH(sdf$RH, es, sdf$p_start),
        R_v = (sdf$delta_vapor/1000 + 1)*ref_ratio
      ) %>% unlist()
      ode(
        y = variables_initial, 
        times = seq(from = 0, to = round(max_height/dz), by = 1), 
        func = calc_derivs, 
        parms = as.list(.)) %>% 
        as_data_frame() 
    }) %>% 
    left_join(parameters, by = "scenario") %>% 
    mutate(alpha = calc_alpha(Temp.K),
           R_c = ifelse(R_v == R_v[1], NA, alpha * R_v)) %>% 
    mutate(`Temp [K]` = Temp.K) %>%
    mutate(`Vapor fraction [%]` = 100*q) %>%
    mutate(`Delta xx vapor [permil]` = 1000*(R_v / ref_ratio - 1),
           `Delta xx condensate [permil]` = 1000*(R_c / ref_ratio - 1)) %>%
    gather(variable, value, contains("[")) %>% 
    mutate(variable = sub("xx", isotope, variable),
           panel = sub("(vapor |condensate )", "", variable)) 
}

# derivatives function
calc_derivs <- function(t, y, parms) {
  with (c(as.list(y), parms), {
    
    # calculate pressure [in Pascal]
    p <- exp(-G) * p_start
    
    # calculate latent heat (Source: wikipedia.org/wiki/Latent_heat)
    Temp.C <- Temp.K - 273.15
    L <- (2500.8 - 2.36 * Temp.C + 0.0016 * Temp.C^2 - 0.00006*Temp.C^3 )*1000 # in J/kg
    
    # calculate saturation vapor pressure and saturatin mass mixing ratio
    des <- calc_es(Temp.K) - es # change in es
    #es <- calc_es(Temp.K) # new es
    qs <- calc_qs(es, p)
    
    # are we condensing?
    condensing <- q >= qs
    
    # scaling (pressure / height scaling)
    dG <- g/(R_air * Temp.K) * dz 
    
    if (condensing) {
      # derivative of the saturation vapor pressure
      ln_es_deriv <- 1/es * calc_des_dT(Temp.K)
      L <- calc_L(Temp.K)
      dT <- -(R_air * Temp.K + L * qs) / (Cp_air + L* qs * ln_es_deriv) * dG 
      dq <- qs - q # drop q down to the saturation (qs)
      dR <- R_v * (calc_alpha(Temp.K) - 1) * 1/q * dq
    } else {
      # not condensing (yay easy)
      dT <- - (R_air/Cp_air) * Temp.K * dG
      dq <- 0
      dR <- 0
    }
    
    # make sure to return in the same order as initial values
    return(list(c(dz, dG, dT, des, dq, dR)))
  })
}