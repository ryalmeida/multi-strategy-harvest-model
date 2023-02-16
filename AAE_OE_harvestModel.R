#'"MULTIPLE CO-OCCURRING BIOECONOMIC DRIVERS OF OVREXPLOITATION ACCELERATE
#'RARE SPECIES EXTINCTION RISK": MODEL CODE
#'
#'AUTHOR: Ryan J. Almeida
#'
#'DATE: 2023/15/02
#'
#'REQUIRED PACKAGES: tidyverse
#'
#'DESCRIPTION: This script contains the source code for the 
#'multi-species harvest model described in Almeida et al. 2023
#'and reproduces the analyses from the paper's numerical example 
#'focused on pangolin harvest.



#Open tidyverse
library(tidyverse)

### MODEL FUNCTIONS ### ----

#'AAE_OE Harvest Model
#'
#'@param state starting values for three state variables (N1, N2, E)
#'@param parms set of parameter values
#'
#'@description ODE-based harvest model that describes the impact of
#'harvest effort (E) on two species populations: N1, a focal species with 
#'a variable price-abundance relationship, and N2, a common species that
#'is sold for a constant p0. Harvest effort (E) in the model is split between
#'two strategies: pursuit harvest, where only the focal species is hunted, and
#'opportunistic harvest, where both the common species and the focal species are
#'hunted. The net profit made from harvesting and selling species in a given time step
#'drives the overall harvest effort in the next time step.
AAE_OE <- function(state, parms){
  
  with(as.list(c(state,parms)), {
    if (N1 >= Ne & N2 >= Ne & E >= 0){ #if biologically feasible, with extinction threshold of Ne
      
      #incorporate variability in harvest strategy (s)
      s <- rnorm(1,s,s.sd) 
      
      #define lower bound of s
      if (s < 0){ 
        s <- 0
      }
      
      #define upper bound of s 
      if (s > 1){ 
        s <- 1
      }
      
      #incorporate variability into p0
      p0 <- rnorm(1,p0,p0.sd) 
      
      
      #ODE describing population dyanmics of focal species
      dN1 <- ((r1*N1*(1-(N1/k1))) - (q1*N1*E*s) - (gamma*N1*E*(1-s))) 
      
      #ODE describing population dynamics of common species
      dN2 <- ((r2*N2*(1-(N2/k2))) - (q2*N2*E*(1-s))) 
      
      #Define price-abundance relationship
      PN <- a + (b/(N1^z + 1))
      
      #ODE describing hunting effort dynamics
      dE <- alpha*((((PN*(q1*N1*E)) - (c*E))*s)
                   + ((1-s)*((PN*gamma*N1*E) + (p0*q2*N2*E) - (c0*E))))
      
      #return state variables
      return((c(dN1, dN2, dE))) 
    }
    else{ #if not biologically feasible, return 0
      return((c(0,0,0)))
    }
  })
}

#'Euler solver for numerical integration of AAE_OE model
#'
#'@param t0 initial time 
#'@param step time step for integration
#'@param N01 initial population size of focal species
#'@param N02 initial population size of common species
#'@param E0 initial harvest effort
#'@param Ne quasi-extinction threshold
#'@param tmax ending time
#'@param params vector of parameter values to pass into AAE_OE model
#'@param g environmental stochasticity constant 
#'
#'@description Euler solver designed to numerically integrate AAE_OE model.
euler_int <- function(t0, step, N01, N02, E0, Ne, tmax, params, g){
  
  #determine time steps to be evaluated
  steps <- seq(t0,tmax, step)
  
  #initialize vector of N1 values
  N1_vals <- rep(0,length(steps))
  #initialize vector of N2 values
  N2_vals <- rep(0,length(steps))
  #set first N1 value to N0
  N1_vals[1] <- N01
  #set first N2 value to N0
  N2_vals[1] <- N02
  #initialize vector of E values
  E_vals <- rep(0,length(steps))
  #set first E value to E0
  E_vals[1] <- E0
  for (i in 2:length(steps)){ #while we haven't reached stopping time
    #set the initial conditions for this time
    state <- c(N1 = N1_vals[i-1], N2 = N2_vals[i-1], E = E_vals[i-1])
    #evaluate dN1, dN2, dE
    fNE <- AAE_OE(state, c(params, Ne = Ne))
    #evaluate N1t with noise
    N1_vals[i] <- N1_vals[i-1] + step*fNE[1] +
      sqrt(step)*g*N1_vals[i-1]*rnorm(1)
    #if N1t is below our threshold
    if (N1_vals[i] < Ne) {
      #we consider it extinct, set to 0
      N1_vals[i] <- 0
      #break loop
      break
    }
    #evaluate Nt with noise
    N2_vals[i] <- N2_vals[i-1] + step*fNE[2] +
      sqrt(step)*g*N2_vals[i-1]*rnorm(1)
    #if Nt is below our threshold
    if (N2_vals[i] < Ne) {
      #we consider it extinct, set to 0
      N2_vals[i] <- 0
      #break loop
      break
    }
    #evaluate Et
    E_vals[i] <- E_vals[i-1] + step*fNE[3]
    #go to next time step
  }
  #return results
  out <- data.frame("t" = steps, "N1" = N1_vals, "N2" = N2_vals, "E" = E_vals)
  return(out)
}

### EXAMPLE SIMULATIONS ### ----

#initialize parameters for ground pangolin case study
pangolins <- c(r1 = 0.1, q1 = 4.2E-5, gamma = 4.2E-6, k1 = 50000)

#initialize parameters for a slow life-history common species
slow.sp <- c(r2 = 0.1, q2 = 4.2E-5, k2 = 100000)

#initialize parameters for a fast life-history common species
fast.sp <- c(r2 = 1, q2 = 4.2E-4, k2 = 500000)

#initialize select bioeconomic parameters
bioecon.parms<- c(alpha = 0.001, p0 = 300, p0.sd = 0, c = 300, c0 = 150)

#initialize constant price scenario parameters
PN.constant <- c(a=300, b = 0, z = 0)

#initialize weak price-abundance relationship parameters
PN.low <- c(a = 50, b = 3.8E4, z = 0.5)

#initialize strong price-abundance relationship parameters
PN.high <- c(a = 50, b = 9.2E8, z = 1.5)

#initialize remaining parameters for integration
tmax <- 5000
step <- 1E-2
N01 <- 24000
N02 <- 50000
Ne <- 10
g <- 0.00
nsim <- 1000

#Slow common species

#Constant price scenario, opportunistic harvest
constant.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                         tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.constant, s = 0, s.sd = 0), g = g)
#Constant price scenario, pursuit harvest
constant.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                         tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.constant, s = 1, s.sd = 0), g = g)
#Constant price scenario, both strategies
constant.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                          tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.constant, s = 0.5, s.sd = 0), g = g)
#Low P-A scenario, opportunistic harvest
low.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                    tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.low, s = 0, s.sd = 0), g = g)
#Low P-A scenario, pursuit harvest
low.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                    tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.low, s = 1, s.sd = 0), g = g)
#Low P-A scenario, both strategies
low.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.low, s = 0.5, s.sd = 0), g = g)
#High P-A scenario, opportunistic harvest
high.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.high, s = 0, s.sd = 0), g = g)
#High P-A scenario, pursuit harvest
high.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.high, s = 1, s.sd = 0), g = g)
#High P-A scenario, both strategies
high.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                      tmax = tmax, params = c(pangolins, slow.sp,  bioecon.parms, PN.high, s = 0.5, s.sd = 0), g = g)


#Reformat output dataframes for plotting
constant.OE.df <- pivot_longer(constant.OE, 2:4, "var")
constant.OE.df$s <- 0
constant.OE.df$PN <- "Constant price"

constant.PH.df <- pivot_longer(constant.PH, 2:4, "var")
constant.PH.df$s <- 1
constant.PH.df$PN <- "Constant price"

constant.mix.df <- pivot_longer(constant.mix, 2:4, "var")
constant.mix.df$s <- 0.5
constant.mix.df$PN <- "Constant price"

low.OE.df <- pivot_longer(low.OE, 2:4, "var")
low.OE.df$s <- 0
low.OE.df$PN <- "Low rarity value"

low.PH.df <- pivot_longer(low.PH, 2:4, "var")
low.PH.df$s <- 1
low.PH.df$PN <- "Low rarity value"

low.mix.df <- pivot_longer(low.mix, 2:4, "var")
low.mix.df$s <- 0.5
low.mix.df$PN <- "Low rarity value"

high.OE.df <- pivot_longer(high.OE, 2:4, "var")
high.OE.df$s <- 0
high.OE.df$PN <- "high rarity value"

high.PH.df <- pivot_longer(high.PH, 2:4, "var")
high.PH.df$s <- 1
high.PH.df$PN <- "high rarity value"

high.mix.df <- pivot_longer(high.mix, 2:4, "var")
high.mix.df$s <- 0.5
high.mix.df$PN <- "high rarity value"

constant.df <- rbind(constant.OE.df,constant.PH.df, constant.mix.df)
constant.df$common <- "Slow"
low.df <- rbind(low.OE.df,low.PH.df, low.mix.df)
low.df$common <- "Slow"
high.df <- rbind(high.OE.df,high.PH.df, high.mix.df)
high.df$common <- "Slow"

#Label harvest strategies according to s value
harvest.strats <- c("Opportunistic harvest", "Pursuit harvest", "Mixed strategies")
names(harvest.strats) <- c(0,1,0.5)

#Constant OH plot
constant.OH.sim <- ggplot(filter(constant.df, s == 0), aes(x = log10(t), y = log10(value+1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.OH.sim


#Constant PH plot

constant.PH.sim <- ggplot(filter(constant.df, s == 1), aes(x = log10(t), y = log10(value +1 ), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.PH.sim

#Constant Mix plot

constant.mix.sim <- ggplot(filter(constant.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.mix.sim

#low OH plot
low.OH.sim <- ggplot(filter(low.df, s == 0), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none")+
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

low.OH.sim



#low PH plot
low.PH.sim <- ggplot(filter(low.df, s == 1), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        # panel.border = element_rect(colour = "black", fill = NA), legend.position = c(0.6,0.8), #uncomment for legend
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 26), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)")

low.PH.sim

#low Mix plot

low.mix.sim <- ggplot(filter(low.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

low.mix.sim

#high OH plot
high.OH.sim <- ggplot(filter(high.df, s == 0), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none")+
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.OH.sim


#high PH plot
high.PH.sim <- ggplot(filter(high.df, s == 1), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        # panel.border = element_rect(colour = "black", fill = NA), legend.position = c(0.6,0.8), #uncomment for legend
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 26), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.PH.sim

#high Mix plot

high.mix.sim <- ggplot(filter(high.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.mix.sim

###############fast common species


constant.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                         tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.constant, s = 0, s.sd = 0), g = g)
constant.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                         tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.constant, s = 1, s.sd = 0), g = g)
constant.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                          tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.constant, s = 0.5, s.sd = 0), g = g)
low.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                    tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.low, s = 0, s.sd = 0), g = g)
low.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                    tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.low, s = 1, s.sd = 0), g = g)
low.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.low, s = 0.5, s.sd = 0), g = g)
high.OE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.high, s = 0, s.sd = 0), g = g)
high.PH <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                     tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.high, s = 1, s.sd = 0), g = g)
high.mix <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                      tmax = tmax, params = c(pangolins, fast.sp,  bioecon.parms, PN.high, s = 0.5, s.sd = 0), g = g)


constant.OE.df <- pivot_longer(constant.OE, 2:4, "var")
constant.OE.df$s <- 0
constant.OE.df$PN <- "Constant price"

constant.PH.df <- pivot_longer(constant.PH, 2:4, "var")
constant.PH.df$s <- 1
constant.PH.df$PN <- "Constant price"

constant.mix.df <- pivot_longer(constant.mix, 2:4, "var")
constant.mix.df$s <- 0.5
constant.mix.df$PN <- "Constant price"

low.OE.df <- pivot_longer(low.OE, 2:4, "var")
low.OE.df$s <- 0
low.OE.df$PN <- "Low rarity value"

low.PH.df <- pivot_longer(low.PH, 2:4, "var")
low.PH.df$s <- 1
low.PH.df$PN <- "Low rarity value"

low.mix.df <- pivot_longer(low.mix, 2:4, "var")
low.mix.df$s <- 0.5
low.mix.df$PN <- "Low rarity value"

high.OE.df <- pivot_longer(high.OE, 2:4, "var")
high.OE.df$s <- 0
high.OE.df$PN <- "high rarity value"

high.PH.df <- pivot_longer(high.PH, 2:4, "var")
high.PH.df$s <- 1
high.PH.df$PN <- "high rarity value"

high.mix.df <- pivot_longer(high.mix, 2:4, "var")
high.mix.df$s <- 0.5
high.mix.df$PN <- "high rarity value"

constant.df <- rbind(constant.OE.df,constant.PH.df, constant.mix.df)
constant.df$common <- "fast"
low.df <- rbind(low.OE.df,low.PH.df, low.mix.df)
low.df$common <- "fast"
high.df <- rbind(high.OE.df,high.PH.df, high.mix.df)
high.df$common <- "fast"

harvest.strats <- c("Opportunistic harvest", "Pursuit harvest", "Mixed strategies")
names(harvest.strats) <- c(0,1,0.5)

#Constant OH plot
constant.OH.sim <- ggplot(filter(constant.df, s == 0), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.OH.sim


#Constant PH plot

constant.PH.sim <- ggplot(filter(constant.df, s == 1), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.PH.sim

#Constant Mix plot

constant.mix.sim <- ggplot(filter(constant.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

constant.mix.sim

#low OH plot
low.OH.sim <- ggplot(filter(low.df, s == 0), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none")+
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

low.OH.sim



#low PH plot
low.PH.sim <- ggplot(filter(low.df, s == 1), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        # panel.border = element_rect(colour = "black", fill = NA), legend.position = c(0.6,0.8), #uncomment for legend
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 26), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)")

low.PH.sim

#low Mix plot

low.mix.sim <- ggplot(filter(low.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

low.mix.sim

#high OH plot
high.OH.sim <- ggplot(filter(high.df, s == 0), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 24, y = 65000, label  = "Opportunistic harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none")+
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.OH.sim


#high PH plot
high.PH.sim <- ggplot(filter(high.df, s == 1), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        # panel.border = element_rect(colour = "black", fill = NA), legend.position = c(0.6,0.8), #uncomment for legend
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), legend.text = element_text(size = 26), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.PH.sim

#high Mix plot

high.mix.sim <- ggplot(filter(high.df, s == 0.5), aes(x = log10(t), y = log10(value + 1), col = var))+
  geom_line(size = 3.5) +
  facet_grid(cols = vars(s), labeller = labeller(s = harvest.strats)) +
  scale_color_manual(labels = c("Harvest effort", "Pangolin population size", "Common species population size"), values = c("orangered2", "navy", "dark green")) +
  #annotate(geom = "text", x = 75, y = 65000, label  = "Pursuit harvest", size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 26, colour = "black"),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.position = "none",
        legend.title = element_blank(), strip.background = element_blank()) +
  ylab("log(Individuals)\n") +
  xlab("\nlog(Years)") 

high.mix.sim

### EQUILIBRIUM ANALYSIS ### ----

#' In this analysis, we determine stationary state values for the focal species
#' and harvest effort across a range of different harvest strategies, under 
#' three different bioeconomic scenarios, and two different common species life history types. 

#Set parameter values for each scenario
pangolins <- c(r1 = 0.1, q1 = 4.2E-5, gamma = 4.2E-6, k1 = 50000)
slow.sp <- c(r2 = 0.1, q2 = 4.2E-5, k2 = 100000)
fast.sp <- c(r2 = 1, q2 = 4.2E-4, k2 = 500000)
bioecon.parms<- c(alpha = 0.001, p0 = 300, p0.sd = 0, c = 300, c0 = 150)
PN.constant <- c(a=300, b = 0, z = 0)
PN.low <- c(a = 50, b = 3.8E4, z = 0.5)
PN.high <- c(a = 50, b = 9.2E8, z = 1.5)
tmax <- 5000
step <- 1E-1
N01 <- 24000
N02 <- 50000
Ne <- 10
g <- 0
nsim <- 1

#set range of s values representing different harvest strategies
s.vals <- seq(0,1,by=0.01)

#initialize data frame for results
equilibrium.df <- data.frame()

for (s in s.vals){ #for each harvest strategy 
  
  #initialize vectors for equilibrium pangolin population size for each bioeconomic and common species scenario
  N1.vals.constant.slow <- rep(0,nsim)
  N1.vals.constant.fast <- rep(0,nsim)
  N1.vals.low.slow <- rep(0,nsim)
  N1.vals.low.fast <- rep(0,nsim)
  N1.vals.high.slow <- rep(0,nsim)
  N1.vals.high.fast<- rep(0,nsim)
  
  #initialize vectors for equilibrium common species population size for each bioeconomic and common species scenario
  N2.vals.constant.slow <- rep(0,nsim)
  N2.vals.constant.fast <- rep(0,nsim)
  N2.vals.low.slow <- rep(0,nsim)
  N2.vals.low.fast <- rep(0,nsim)
  N2.vals.high.slow <- rep(0,nsim)
  N2.vals.high.fast<- rep(0,nsim)
  
  #initialize vectors for equilibrium harvest effort for each bioeconomic and common species scenario
  E.vals.constant.slow <- rep(0,nsim)
  E.vals.low.slow <- rep(0,nsim)
  E.vals.high.slow <- rep(0,nsim)
  E.vals.constant.fast <- rep(0,nsim)
  E.vals.low.fast <- rep(0,nsim)
  E.vals.high.fast <- rep(0,nsim)
  
  for (i in 1:nsim){
    
    #run simulations for each bioeconomic & life-history scenario
    constant.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                                   tmax = tmax, params = c(pangolins, bioecon.parms, slow.sp, PN.constant, s = s, s.sd = 0), g = g)
    constant.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                                   tmax = tmax, params = c(pangolins, bioecon.parms, fast.sp, PN.constant, s = s, s.sd = 0), g = g)
    low.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne,
                              tmax = tmax, params = c(pangolins, bioecon.parms, slow.sp, PN.low, s = s, s.sd = 0), g = g)
    low.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne,
                              tmax = tmax, params = c(pangolins, bioecon.parms, fast.sp, PN.low, s = s, s.sd = 0), g = g)
    high.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                               tmax = tmax, params = c(pangolins, bioecon.parms, slow.sp, PN.high, s = s, s.sd = 0), g = g)
    high.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                               tmax = tmax, params = c(pangolins, bioecon.parms, fast.sp, PN.high, s = s, s.sd = 0), g = g)
    
    #for numerical simulations, we will define equilibrium values as average from last 10% of time steps
    t.start <- (0.9*tmax/step) + 1
    t.end <- tmax/step+1
    
    #Store each variable's equilibrium state when common species has slow life history
    N1.vals.constant.slow[i] <- mean(constant.slow.out$N1[t.start:t.end])
    E.vals.constant.slow[i] <- mean(constant.slow.out$E[t.start:t.end])
    N1.vals.low.slow[i] <- mean(low.slow.out$N1[t.start:t.end])
    E.vals.low.slow[i] <- mean(low.slow.out$E[t.start:t.end])
    N1.vals.high.slow[i] <- mean(high.slow.out$N1[t.start:t.end])
    E.vals.high.slow[i] <- mean(high.slow.out$E[t.start:t.end])
    N2.vals.constant.slow[i] <- mean(constant.slow.out$N2[t.start:t.end])
    N2.vals.low.slow[i] <- mean(low.slow.out$N2[t.start:t.end])
    N2.vals.high.slow[i] <- mean(high.slow.out$N2[t.start:t.end])
    
    #Store each variable's equilibrium state when common species has fast life history
    N1.vals.constant.fast[i] <- mean(constant.fast.out$N1[t.start:t.end])
    E.vals.constant.fast[i] <- mean(constant.fast.out$E[t.start:t.end])
    N1.vals.low.fast[i] <- mean(low.fast.out$N1[t.start:t.end])
    E.vals.low.fast[i] <- mean(low.fast.out$E[t.start:t.end])
    N1.vals.high.fast[i] <- mean(high.fast.out$N1[t.start:t.end])
    E.vals.high.fast[i] <- mean(high.fast.out$E[t.start:t.end])
    N2.vals.constant.fast[i] <- mean(constant.fast.out$N2[t.start:t.end])
    N2.vals.low.fast[i] <- mean(low.fast.out$N2[t.start:t.end])
    N2.vals.high.fast[i] <- mean(high.fast.out$N2[t.start:t.end])
    
  }
  #store results into output dataframe
  constant.slow.df <- data.frame("N1" = mean(N1.vals.constant.slow), "N2" = mean(N2.vals.constant.slow), "E" = mean(E.vals.constant.slow), 
                                 "s" = s, "PN" = "Constant price", "common" = "Slow")
  low.slow.df <- data.frame("N1" = mean(N1.vals.low.slow), "N2" = mean(N2.vals.low.slow), "E" = mean(E.vals.low.slow), 
                            "s" = s, "PN" = "Low rarity value", "common" = "Slow")
  high.slow.df <- data.frame("N1" = mean(N1.vals.high.slow), "N2" = mean(N2.vals.high.slow), "E" = mean(E.vals.high.slow), 
                             "s" = s, "PN" = "High rarity value", 'common' = "Slow")
  constant.fast.df <- data.frame("N1" = mean(N1.vals.constant.fast), "N2" = mean(N2.vals.constant.fast), "E" = mean(E.vals.constant.fast), 
                                 "s" = s, "PN" = "Constant price", "common" = "Fast")
  low.fast.df <- data.frame("N1" = mean(N1.vals.low.fast), "N2" = mean(N2.vals.low.fast), "E" = mean(E.vals.low.fast), 
                            "s" = s, "PN" = "Low rarity value", "common" = "Fast")
  high.fast.df <- data.frame("N1" = mean(N1.vals.high.fast), "N2" = mean(N2.vals.high.fast), "E" = mean(E.vals.high.fast), 
                             "s" = s, "PN" = "High rarity value", 'common'= "Fast")
  equilibrium.df <- rbind(equilibrium.df, constant.slow.df, low.slow.df,high.slow.df, constant.fast.df, low.fast.df,high.fast.df)
}

#reformat output data for plotting
equilibrium.df.long <- pivot_longer(equilibrium.df, 1:2, "species")

#plot Figure 2
eq2sp.plot <- ggplot(equilibrium.df.long, aes(x = E, y = log10(value), shape = species, color = s)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  geom_point(size = 4) +
  facet_grid(common ~ PN) +
  scale_color_gradient(low = "navy", high = "orange") +
  scale_shape_manual(values = c(8,19)) +
  guides(shape = "none") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_blank(), 
        legend.position = c(0.05, 0.8),
        legend.title = element_blank()) +
  xlab("Harvest effort") +
  ylab("Population size") 

#View Figure 2
eq2sp.plot


### EXTINCTION PROBABILITY DYNAMICS ### -----

#Set parameter values for analysis
pangolins <- c(r1 = 0.1, q1 = 4.2E-5, gamma = 4.2E-6, k1 = 50000)
slow.sp <- c(r2 = 0.1, q2 = 4.2E-5, k2 = 50000)
fast.sp <- c(r2 = 1, q2 = 4.2E-4, k2 = 500000)
bioecon.parms<- c(alpha = 0.001, p0 = 300, p0.sd = 50, c = 300, c0 = 150)
PN.constant <- c(a=300, b = 0, z = 0)
PN.low <- c(a = 50, b = 3.8E4, z = 0.5)
PN.high <- c(a = 50, b = 9.2E8, z = 1.5)
#Set short time scale of 50 years
tmax <- 50
step <- 1E-1
N01 <- 24000
N02 <- 40000
Ne <- 1
g <- 0.05
s.sd <- 0.15 
#set number of simulations
nsim <- 1000

#Initialize vector of s values representing different harvest strategies
s.vals <- seq(0,1,by = 0.01)
#Initialize output dataframe
Pe.vals <- data.frame()

for (s in s.vals){ #for each harvest strategy represented by s
  #Print current s value to console to track progress
  print(paste(c("Current s:", s)))
  
  #Initialize vectors to store variable outputs for different scenarios
  N1.vals.slow.low <- rep(N01,nsim)
  t.vals.slow.low <- rep(NA,nsim)
  N1.vals.slow.high <- rep(N01,nsim)
  t.vals.slow.high <- rep(NA,nsim)
  N1.vals.slow.constant <- rep(N01,nsim)
  t.vals.slow.constant <- rep(NA,nsim)
  N1.vals.fast.low <- rep(N01,nsim)
  t.vals.fast.low <- rep(NA,nsim)
  N1.vals.fast.high <- rep(N01,nsim)
  t.vals.fast.high <- rep(NA,nsim)
  N1.vals.fast.constant <- rep(N01,nsim)
  t.vals.fast.constant <- rep(NA,nsim)
  
  for (i in 1:nsim){ #for each simulation
    
    #Simulate 50 years of harvest for weak price-abundance relationship
    low.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                              tmax = tmax, params = c(pangolins, PN.low, bioecon.parms, slow.sp, s = s, s.sd = s.sd), g = g) 
    low.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                              tmax = tmax, params = c(pangolins, PN.low, bioecon.parms, fast.sp, s = s, s.sd = s.sd), g = g) 
    
    #Store population size at the end of the simulation
    current.N.slow <- low.slow.out$N1[length(low.slow.out$N1)]
    current.N.fast <- low.fast.out$N1[length(low.fast.out$N1)]
    
    #determine which time step the species went extinct
    time.extinct.low.slow <- which(low.slow.out$N1 == 0)
    
    if(length(time.extinct.low.slow) > 0){ #if the species was driven extinct
      t.vals.slow.low[i] <- low.slow.out$t[time.extinct.low.slow[1]] #store time of extinction
    }
    
    if (current.N.slow < N01){ #if the species' population was depleted
      N1.vals.slow.low[i] <- current.N.slow #determine final population size
    }
    
    #Follow procedure above but for constant price scenario
    constant.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                                   tmax = tmax, params = c(pangolins, bioecon.parms, slow.sp, PN.constant, s = s, s.sd), g = g)
    constant.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                                   tmax = tmax, params = c(pangolins, bioecon.parms, fast.sp, PN.constant, s = s, s.sd), g = g) 
    
    current.N.slow <- constant.slow.out$N1[length(constant.slow.out$N1)]
    current.N.fast <- constant.fast.out$N1[length(constant.fast.out$N1)]
    
    time.extinct.constant.slow <- which(constant.slow.out$N1 == 0)
    if(length(time.extinct.constant.slow) > 0){
      t.vals.slow.constant[i] <- constant.slow.out$t[time.extinct.constant.slow[1]]
    }
    
    if (current.N.slow < N01){
      N1.vals.slow.constant[i] <- current.N.slow
    }
    
    #Follow procedure above but for strong price-abundance relationship scenario
    high.slow.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                               tmax = tmax, params = c(pangolins, bioecon.parms, slow.sp, PN.high, s = s, s.sd), g = g)
    high.fast.out <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, Ne = Ne, 
                               tmax = tmax, params = c(pangolins, bioecon.parms, fast.sp, PN.high, s = s, s.sd), g = g) 
    
    current.N.slow <- high.slow.out$N1[length(high.slow.out$N1)]
    current.N.fast <- high.fast.out$N1[length(high.fast.out$N1)]
    
    time.extinct.high.slow <- which(high.slow.out$N1 == 0)
    if(length(time.extinct.high.slow) > 0){
      t.vals.slow.high[i] <- high.slow.out$t[time.extinct.high.slow[1]]
    }
    
    if (current.N.slow < N01){
      N1.vals.slow.high[i] <- current.N.slow
    }
  }
  
  #Store results in output dataframe
  slow.constant.df <- data.frame("PN" = "Constant price", "depletion" = (mean(N1.vals.slow.constant/N01)), 
                                 Pe = sum(N1.vals.slow.constant == 0)/nsim, "t" = mean(t.vals.slow.constant), "common" = "Slow")
  slow.low.df <- data.frame("PN" = "Low rarity value", "depletion" = (mean(N1.vals.slow.low/N01)), Pe = sum(N1.vals.slow.low == 0)/nsim,
                            "t" = mean(t.vals.slow.low), "common" = "Slow")
  slow.high.df <- data.frame("PN" = "High rarity value", "depletion" = (mean(N1.vals.slow.high/N01)), Pe = sum(N1.vals.slow.high == 0)/nsim,
                             "t" = mean(t.vals.slow.high), "common" = "Slow")
  fast.constant.df <- data.frame("PN" = "Constant price", "depletion" = (mean(N1.vals.fast.constant/N01)), 
                                 Pe = sum(N1.vals.fast.constant == 0)/nsim, "t" = mean(t.vals.fast.constant), "common" = "Fast")
  fast.low.df <- data.frame("PN" = "Low rarity value", "depletion" = (mean(N1.vals.fast.low/N01)), Pe = sum(N1.vals.fast.low == 0)/nsim,
                            "t" = mean(t.vals.fast.low), "common" = "Fast")
  fast.high.df <- data.frame("PN" = "High rarity value", "depletion" = (mean(N1.vals.fast.high/N01)), Pe = sum(N1.vals.fast.high == 0)/nsim,
                             "t" = mean(t.vals.fast.high), "common" = "Fast")
  current <- rbind(slow.constant.df, slow.low.df, slow.high.df, fast.constant.df, fast.low.df, fast.high.df)
  current$s <- s
  Pe.vals <- rbind(Pe.vals, current)
}

#Reformat dataframe for plotting
Pe.vals.long <- pivot_longer(Pe.vals, 2:3,names_to = "metric")

#Define colorblind friendly scale
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Plot figure 3
PeFig <- ggplot(Pe.vals.long, aes(x = s, y = value, colour = metric, linetype = common)) +
  #geom_jitter(size = 3) +
  geom_line(size = 3) +
  facet_grid(cols = vars(PN)) +
  guides(linetype = "none") +
  scale_color_manual(values = colorBlindBlack8, labels = c("Proportion population\nremaining", "Extinction probability")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), legend.title = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_blank(), 
        legend.position = c(0.2,0.2),
        legend.text = element_text(size = 12, colour = "black")) +
  ylab("") +
  xlab("\nProportion pursuit hunting (s)")

#View figure 3
PeFig

### EXTINCTION PROBABILITY SENSITIVITY ANALYSIS ### -----

#' In this analysis, we determine how sensitive focal species extinction probability (Pe)
#' is to major changes in parameter values. For each baseline bioeconomic parameter value, 
#' we multiply that parameter by scalars ranging from 0.001 to 1000 and determine how Pe
#' changes. 

#set initial conditions
tmax <- 50
step <- 1E-1
N01 <-24000
N02 <- 500000
g <- 0.05
Ne <- 1
nsim <- 10
s_vals <- seq(0,1,by = 0.01)

#initialize data frame for results
parms.df <- data.frame()

#define parameter names
parm_names <- c("gamma", "alpha", "b", "z", "k1", "c", "c0", "a", "r1", "p0", "q1", "p0.sd", "Ne", "g", "r2", "k2", "q2")
#define scalar multipliers to alter parameter values
parm_range <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
#set initial parameter values 
parms_base <- c(r1 = 0.1, q1 = 4.2E-5, gamma = 4.2E-6, k1 = 50000, alpha = 0.01, a = 50, p0 = 300, p0.sd = 50, c = 300, c0 = 150, b = 3.8E4, z = 0.5, s.sd = 0.15, r2 = 0.1, k2 = 50000, q2 = 4.2E-5)



for (i in 1:length(parm_names)){ #for each parameter
  current_parm <- parm_names[i] #store the current parameter
  pos <- which(names(parms_base) == current_parm) #determine which position the parameter is in
  if (current_parm == "g"){
    for (j in 1:length(parm_range)){ #for each new value of parameter
      g.current <- parm_range[j]*g
      print(paste("Current parameter:",current_parm, ", current val:",g.current)) #print which parameter value we're on
      for (s in s_vals){ #for each s
        N1.vals <- rep(0,nsim) #initialize vector to store pangolin pop size
        E_vals <- rep(0,nsim) #initialize vector to store harvest efort
        for (k in 1:nsim){ #for each simulation
          out.AAE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, #run simulation
                               tmax = tmax, params = parms_base, g = g.current, Ne = Ne)
          N1.vals[k] <- out.AAE$N1[length(out.AAE$N1)] #save N val
          E_vals[k] <- out.AAE$E[length(out.AAE$E)] #save E val
        }
        df.N <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(N1.vals), "var" = "N", "Pe" = (sum(N1.vals == 0)/nsim))
        df.E <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(E_vals),  "var" = "E", "Pe" = (sum(E_vals == 0)/nsim))
        parms.df <- rbind(parms.df,df.N, df.E)
      }
    }
    next
  }
  if (current_parm == "Ne"){
    for (j in 1:length(parm_range)){ #for each new value of parameter
      Ne.current <- parm_range[j]*Ne
      print(paste("Current parameter:",current_parm, ", current val:",Ne.current)) #print which parameter value we're on
      for (s in s_vals){ #for each s
        N1.vals <- rep(0,nsim) #initialize vector to store pangolin pop size
        E_vals <- rep(0,nsim) #initialize vector to store harvest efort
        for (k in 1:nsim){ #for each simulation
          out.AAE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, #run simulation
                               tmax = tmax, params = parms_base, g = g, Ne = Ne.current)
          N1.vals[k] <- out.AAE$N1[length(out.AAE$N1)] #save N val
          E_vals[k] <- out.AAE$E[length(out.AAE$E)] #save E val
        }
        df.N <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(N1.vals), "var" = "N", "Pe" = (sum(N1.vals == 0)/nsim))
        df.E <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(E_vals),  "var" = "E", "Pe" = (sum(E_vals == 0)/nsim))
        parms.df <- rbind(parms.df,df.N, df.E)
      }
    }
    next
  }
  for (j in 1:length(parm_range)){ #for each new value of parameter
    current_set <- parms_base #store the original parameter values
    current_set[pos] <- parms_base[pos]*parm_range[j] #multiply the parameter by the scalar
    print(paste("Current parameter:",current_parm, ", current val:",current_set[pos])) #print which parameter value we're on
    for (s in s_vals){ #for each s
      N1.vals <- rep(0,nsim) #initialize vector to store pangolin pop size
      E_vals <- rep(0,nsim) #initialize vector to store harvest efort
      for (k in 1:nsim){ #for each simulation
        out.AAE <- euler_int(t0 = 0, step = step, N01 = N01, N02 = N02, E0 = 1, #run simulation
                             tmax = tmax, params = current_set, Ne = Ne, g = g)
        N1.vals[k] <- out.AAE$N1[length(out.AAE$N1)] #save N val
        E_vals[k] <- out.AAE$E[length(out.AAE$E)] #save E val
      }
      df.N <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(N1.vals), "var" = "N", "Pe" = (sum(N1.vals == 0)/nsim))
      df.E <- data.frame("s" = s, "parm_val" = parm_range[j], "parm_name" = current_parm, "val" = mean(E_vals),  "var" = "E", "Pe" = (sum(E_vals == 0)/nsim))
      parms.df <- rbind(parms.df,df.N, df.E)
    }
  }
}

#Re-label parameter levels for clarity when plotting
parms.df$parm_name <- factor(parms.df$parm_name, 
                             levels = c("r1", "k1", "c0", "gamma", "p0", "p0.sd", 
                                        "q1", "c", "a", "b", "z", "alpha", "Ne", "g", 
                                        "r2", "k2", "q2"))

#Plot sensitivity analysis results without Ne or g parameters
parm.Pe.fig <- ggplot(filter(parms.df, parm_name != "Ne", parm_name != "g"), aes(x = s, y = as.factor(parm_val), fill = Pe)) +
  geom_tile() +
  facet_wrap(~parm_name,
             labeller = labeller(parm_name = c("r1" = "Pangolin growth rate (r1)","a" = "Minimum selling price (a)", "alpha" = "Harvest reactivity to\n market dynamics (alpha)", 
                                               "b" = "Selling price\n at extinction (b)", "c" = "Pursuit harvest cost (c)" , "c0" = "Opportunistic\n harvest cost (c0)", "gamma" = "Opportunistic harvest\n catchability (gamma)",
                                               "k1" = "Carrying capacity (k1)", "p0" = "Income from\n other species (p0)", "p0.sd" = "Variation in alternative\n income (p0.sd)", "q1" = "Pursuit harvest\n catchability (q1)",
                                               "z" = "Price-abundance\n scaling parameter (z)", 
                                               "r2" = "Common species\n growth rate (r2)", "k2" = "Common species\n carrying capacity (k2)", "q2" = "Common species\n catchability (q2)")))+
  scale_fill_continuous(low = "white", high = 'black', name = "Extinction\nprobability\n") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        axis.title = element_text(color = "black", size = 36),
        axis.text = element_text(size = 28, colour ="black"),
        strip.text = element_text(color = "black",size = 30),
        legend.title = element_text(color = "black",size = 32),
        legend.text = element_text(color = "black", size = 28),
        panel.spacing = unit(3, "lines"),
        legend.key.size = unit(2, "cm")) +
  ylab("Ratio of new parameter value\n to baseline value") +
  xlab("\nProportion pursuit harvest (s)") 

#View sensitivity analysis results
parm.Pe.fig


