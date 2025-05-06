
#Load packages
source("multilayerHyperNets_simFunctions.R")
library(igraph)
library(data.table)
library(foreach)
library(doParallel)
library(stringr)

#Set number of cores for parallelization
registerDoParallel(cores = 20)

#Create folders in which to store simulation results
run_ID=strftime(Sys.time(), format="d3%Y%m%d%H%M%S")
sim_indData="Sim_individual_level_data"
sim_summaryData="Sim_summary_data"
sim_diffSummaryData = "Sim_diffusion_summary_data"
sim_netData_dyadic = "Sim_dyadic_network"
sim_netData_HO1 = "Sim_higherOrder_network_1"
sim_netData_HO2 = "Sim_higherOrder_network_2"
sim_netData_HO3 = "Sim_higherOrder_network_3"

if(!file.exists(sim_indData)) dir.create(sim_indData)
if(!file.exists(sim_summaryData)) dir.create(sim_summaryData)
if(!file.exists(sim_diffSummaryData)) dir.create(sim_diffSummaryData)
if(!file.exists(sim_netData_dyadic)) dir.create(sim_netData_dyadic)
if(!file.exists(sim_netData_HO1)) dir.create(sim_netData_HO1)
if(!file.exists(sim_netData_HO2)) dir.create(sim_netData_HO2)
if(!file.exists(sim_netData_HO3)) dir.create(sim_netData_HO3)

##Set seed for reproducibility
set.seed(01092025)

#Set social transmission rate; determines baseline probability of learning from an active neighbor in the dyadic layer
social_trans <- 0.1

#Possible radii for determining subgroups
rad <- c(0.15, 0.25, 0.35)

#Possible functions for determining response to dominant individuals
#Only include either suppressive or motivating functions for a run and comment the others out

#Suppressive functions:
domResp <- c(domResponse_Linear_HO, domResponse_Sigmoid_HO, domResponse_Linear_dyadic, domResponse_Sigmoid_dyadic, 
             domResponse_wLinear_HO, domResponse_wSigmoid_HO, domResponse_wLinear_dyadic, domResponse_wSigmoid_dyadic)

#Motivating functions:
#domResp <- c(domResponse_Linear_HO_Increasing, domResponse_Sigmoid_HO_Increasing, domResponse_Linear_dyadic_Increasing, domResponse_Sigmoid_dyadic_Increasing, 
#             domResponse_wLinear_HO_Increasing, domResponse_wSigmoid_HO_Increasing, domResponse_wLinear_dyadic_Increasing, domResponse_wSigmoid_dyadic_Increasing)

#Set name for performance function variants
domRespNames <- c("Linear_HO", "Sigmoid_HO", "Linear_dyadic", "Sigmoid_dyadic",
                  "wLinear_HO", "wSigmoid_HO", "wLinear_dyadic", "wSigmoid_dyadic")

#Distributions for dominance ranks
domDist <- c("domUni", "domExp")

###################################################

#Set number of sims per condition
nSims = 250

foreach(s = 1:nSims) %dopar% {
  
  ##Generate individuals in family/kinship groups
  ind_data <- generate_population(n_families = 25, meanFamilySize = 4, clustering = 0.075, nInitInformed = 1, clusterByFamily = FALSE)
  n_indivs <- nrow(ind_data)
  ind_data$simID <- s
  
  #Set dominance values
  #Exponential distribution is normalized
  #Dominance values are assigned to maintain relative differences in dominance ranks across distributions (e.g., the same individual is assigned the highest value for both distributions, etc.)
  ind_data$domUni <- runif(nrow(ind_data), min = 0, max = 1)
  domExp <- sort(rexp(nrow(ind_data), rate = 2))
  domExp <- sapply(rank(ind_data$domUni), function(x) domExp[x])
  ind_data$domExp <- (domExp - min(domExp))/(max(domExp) - min(domExp))
  ind_data$domExp <- ifelse(ind_data$domExp == 0, min(ind_data[which(ind_data$domExp > 0),]$domExp - 0.0001), ind_data$domExp)
  
  #Keep copy of original population data to use across simulation runs
  ind_dataOrig <- ind_data
  
  ##Create spatial proximity network and hypergraphs
  netList <- generate_latent_space_multilayer_hypergraph(ind_data = ind_data, r = rad)
  
  #Record dyadic- and hypernetworks used for the current set of conditions
  write.csv(netList[[1]], file.path(sim_netData_dyadic, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[2]], file.path(sim_netData_HO1, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[3]], file.path(sim_netData_HO2, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  write.csv(netList[[4]], file.path(sim_netData_HO3, sprintf("simData_%s_%.01i.csv", run_ID, s)))
  
  #Begin running through the different conditions (e.g., dominance distribution, performance functions (i.e., domResp))
  for(d in 1:length(domDist)) {
    
    for(p in 1:length(domResp)) {
      
      #Procedure differs slightly depending on whether performance decisions are made in the dyadic or hypernetwork layers
      if(str_split(domRespNames[p], "_")[[1]][2] == "dyadic") {
        
        #Reset individual-level data to initial values
        ind_data <- ind_dataOrig
        
        #Set dominance distribution and performance function
        domDistribution <- domDist[d]
        
        ind_data$domResponse <- domRespNames[p]
        ind_data$distUsed <- domDistribution
        ind_data$groupRadius <- 0
        
        focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
        domResponseFunction <- domResp[[p]]
        radiusForDomFunction <- 1
        
        ind_data$firstProd <- 0
        
        #Create list for holding output on each time step
        dataList <- vector("list", 8000)
        
        #Set initial time
        t = 1
        
        repeat{
          
          #Identify current set of informed nodes
          informedNodes <- ind_data[which(ind_data$informed == 1),]$id
          
          #Each informed individual determines its likelihood of producing a novel behavior:
          produceTemp <- domResponseFunction(ind_data = ind_data, netList = netList, radius = radiusForDomFunction, informedNodes = informedNodes, domValues = focaldomDist)
          ind_data$produce <- produceTemp
          
          #Determine which informed nodes are producing the novel behavior this time step
          probVect <- runif(nrow(ind_data), min = 0, max = 1)
          ind_data$active <- 0
          for(i in informedNodes) {
            ind_data$active[i] <- ifelse(ind_data$produce[i] >  probVect[i], 1, 0)
          }
          
          #Record first production time for any active demonstrators that have not performed the behavior previously
          for(i in 1:nrow(ind_data)) {
            if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
              ind_data$firstProd[i] <- t
            }
          }
          
          #Determine individuals' likelihood of learning the trait based on their relative connection strength to informed individuals in the dyadic layer
          ind_data$acqProb <- sapply(1:nrow(ind_data), function(x) social_trans * sum(netList[[1]][,x] * ind_data$informed) / sum(netList[[1]][,x]))
          ind_data$acqProb <- ifelse(is.na(ind_data$acqProb), 0, ind_data$acqProb)
          
          #Record those individuals that acquired the trait and note the time step at which this occurred
          acqVect <- runif(nrow(ind_data), min = 0, max = 1)
          ind_data$learned <- 0
          ind_data$learned <- ifelse(ind_data$acqProb > acqVect, 1, 0)
          ind_data$acqTime <- ifelse(ind_data$learned == 1 & ind_data$informed == 0, t, ind_data$acqTime)
          
          #Add newly informed individuals to the list of informed individuals
          ind_data$informed <- ifelse(ind_data$informed == 0, 
                                      ifelse(ind_data$learned == 1, 1, 0), 1)
          
          #Record summary data for the current time step
          dataTemp <- data.table("simID" = s,
                                 "domDist" = domDistribution,
                                 "domResponse" = domRespNames[p],
                                 "groupRadius" = rad[r],
                                 "nInformed" = 0,
                                 "percInformed" = 0,
                                 "nActive" = 0,
                                 "percActive" = 0
          )
          
          dataTemp$nInformed <- sum(ind_data$informed)
          dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
          
          dataTemp$nActive <- sum(ind_data$active)
          dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
          
          dataList[[t]] <- dataTemp
          
          #If all individuals are informed or an inordinately long time has elapsed, end the sim
          if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 8000) {
            break
          }
          
          #If the simulation is continuing, advance to the next time step
          t <- t + 1
        }
        
        #Created aggregated record of diffusion data across time steps
        dataCombined <- rbindlist(dataList)
        dataCombined$timeStep <- 1:nrow(dataCombined)
        
        #Create dataframe for diffusion outcomes
        diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                       "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                       "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                       "numIndivs" = n_indivs, "numProducers" = 0,
                                       "orderDiv" = 0,
                                       "propDiv" = 0,
                                       "timeDelay" = 0)
        
        #Number of individuals that produced the behavior at least once
        numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
        #The mean divergence in rank orders of learning vs. first performance
        orderDiv <- (1/(numProducers - 1)) * 
          sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                    rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
        #The proportion of individuals that differ in their rank orders of learning vs first performance
        propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                         rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
        #The mean delay in time steps between learning the behavior and first performing it
        timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                  ind_data[which(ind_data$firstProd > 0),]$firstProd))
        diffusionSummary$numProducers <- numProducers
        diffusionSummary$orderDiv <- orderDiv
        diffusionSummary$propDiv <- propDiv
        diffusionSummary$timeDelay <- timeDelay
        
        #Output diffusion data
        fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
        
      } else{
        
        #Procedure for higher-order network performance function variants
        if(str_split(domRespNames[p], "_")[[1]][2] == "HO") {
          for(r in 1:length(rad)) {
            
            #Reset individual-level data to initial values
            ind_data <- ind_dataOrig
            
            #Set dominance distribution and performance function
            domDistribution <- domDist[d]
            
            ind_data$domResponse <- domRespNames[p]
            ind_data$distUsed <- domDistribution
            ind_data$groupRadius <- rad[r]
            
            focaldomDist <- as.vector(unlist(ind_data[domDistribution]))
            domResponseFunction <- domResp[[p]]
            radiusForDomFunction <- r + 1
            
            ind_data$firstProd <- 0
            
            #Create list for holding output on each time step
            dataList <- vector("list", 8000)
            
            #Set initial time
            t = 1
            
            repeat{
              
              #Identify current set of informed nodes
              informedNodes <- ind_data[which(ind_data$informed == 1),]$id
              
              #Each informed individual determines its likelihood of producing a novel behavior:
              produceTemp <- domResponseFunction(ind_data = ind_data, netList = netList, radius = radiusForDomFunction, informedNodes = informedNodes, domValues = focaldomDist)
              ind_data$produce <- produceTemp
              
              #Determine which informed nodes are producing the novel behavior this time step
              probVect <- runif(nrow(ind_data), min = 0, max = 1)
              ind_data$active <- 0
              
              for(i in informedNodes) {
                ind_data$active[i] <- ifelse(ind_data$produce[i] >  probVect[i], 1, 0)
              }
              
              #Record first production time for any active demonstrators that have not performed the behavior previously
              for(i in 1:nrow(ind_data)) {
                if(ind_data$active[i] == 1 & ind_data$firstProd[i] == 0) {
                  ind_data$firstProd[i] <- t
                }
              }
              
              #Determine individuals' likelihood of learning the trait based on their relative connection strength to informed individuals in the dyadic layer
              ind_data$acqProb <- sapply(1:nrow(ind_data), function(x) social_trans * sum(netList[[1]][,x] * ind_data$informed) / sum(netList[[1]][,x]))
              ind_data$acqProb <- ifelse(is.na(ind_data$acqProb), 0, ind_data$acqProb)
              
              #Record those individuals that acquired the trait and note the time step at which this occurred
              acqVect <- runif(nrow(ind_data), min = 0, max = 1)
              ind_data$learned <- 0
              ind_data$learned <- ifelse(ind_data$acqProb > acqVect, 1, 0)
              ind_data$acqTime <- ifelse(ind_data$learned == 1 & ind_data$informed == 0, t, ind_data$acqTime)
              
              #Add newly informed individuals to the list of informed individuals
              ind_data$informed <- ifelse(ind_data$informed == 0, 
                                          ifelse(ind_data$learned == 1, 1, 0), 1)
              
              #Record summary data for the current time step
              dataTemp <- data.table("simID" = s,
                                     "domDist" = domDistribution,
                                     "domResponse" = domRespNames[p],
                                     "groupRadius" = rad[r],
                                     "nInformed" = 0,
                                     "percInformed" = 0,
                                     "nActive" = 0,
                                     "percActive" = 0)
              
              dataTemp$nInformed <- sum(ind_data$informed)
              dataTemp$percInformed <- dataTemp$nInformed/nrow(ind_data)
              
              dataTemp$nActive <- sum(ind_data$active)
              dataTemp$percActive <- dataTemp$nActive/dataTemp$nInformed
              
              dataList[[t]] <- dataTemp
              
              #If all individuals are informed or an inordinately long time has elapsed, end the sim
              if(sum(ind_data$firstProd >0) == nrow(ind_data) | t == 8000) {
                break
              }
              
              #If the simulation is continuing, advance to the next time step
              t <- t + 1
            }
            
            #Created aggregated record of diffusion data across time steps
            dataCombined <- rbindlist(dataList)
            dataCombined$timeStep <- 1:nrow(dataCombined)
            
            #Create dataframe for diffusion outcomes
            diffusionSummary <- data.table("simID" = s, "domDist" = domDistribution,
                                           "domResponse" = domRespNames[p], "groupRadius" = rad[r],
                                           "TTD" = max(ind_data$acqTime), "TTFP" = max(ind_data$firstProd),
                                           "numIndivs" = n_indivs, "numProducers" = 0,
                                           "orderDiv" = 0,
                                           "propDiv" = 0,
                                           "timeDelay" = 0)
            
            #Number of individuals that produced the behavior at least once
            numProducers <- nrow(ind_data[which(ind_data$firstProd > 0),])
            #The mean divergence in rank orders of learning vs. first performance
            orderDiv <- (1/(numProducers - 1)) * 
              sum(abs(rank(ind_data[which(ind_data$acqTime > 0  & ind_data$firstProd > 0),]$acqTime) - 
                        rank(ind_data[which(ind_data$acqTime > 0 & ind_data$firstProd > 0),]$firstProd)))
            #The proportion of individuals that differ in their rank orders of learning vs first performance
            propDiv <- sum(rank(ind_data[which(ind_data$firstProd > 0),]$acqTime) != 
                             rank(ind_data[which(ind_data$firstProd > 0),]$firstProd)) / numProducers
            #The mean delay in time steps between learning the behavior and first performing it
            timeDelay <- (1/numProducers) * sum(abs(ind_data[which(ind_data$firstProd > 0),]$acqTime - 
                                                      ind_data[which(ind_data$firstProd > 0),]$firstProd))
            
            diffusionSummary$numProducers <- numProducers
            diffusionSummary$orderDiv <- orderDiv
            diffusionSummary$propDiv <- propDiv
            diffusionSummary$timeDelay <- timeDelay
            
            #Output diffusion data
            fwrite(ind_data, file = file.path(sim_indData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
            fwrite(dataCombined, file = file.path(sim_summaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
            fwrite(diffusionSummary, file = file.path(sim_diffSummaryData, sprintf("simData_%s_%.01i_%.03f_%s_%s.csv", run_ID, s, rad[r], domRespNames[p], domDistribution)))
          }
        }
      }
    }
  }
}

###############################################

