library(tidyverse)
library(lubridate)
library(bnlearn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(networkD3)
library(Hmisc)
library(MASS)
library(visNetwork)




plot.network <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arc_str[,1], structure$arc_str[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  edges <- data.frame(from = structure$arc_str[,1],
                      to = structure$arc_str[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "red",
                      width = abs(structure$arc_str[,3])/10)
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}


set.seed(1020)


data.all <- read.csv2("causal_data1.csv")









data.ext <- data.frame(cases = data.all$cumulative_cases,
                       vacc = data.all$vacc_fraq,
                       HCA_supply = data.all$HCA_supply,
                       HCA_demand = data.all$HCA_demand,
                       incidence = data.all$incidence,
                       death = ifelse(is.na(data.all$deaths_), 0, data.all$deaths_),
                       population_size = data.all$population_size,
                       poulation_density = 1.0*data.all$poulation_density,
                       income = 1.0*data.all$income)

input_vars <- c("cases", "vacc", "HCA_supply", "HCA_demand")
input.vars <- c("cases", "vacc", "HCA_supply", "HCA_demand")



data.comb.iter <- do.call(cbind, combn(input_vars, 2, 
                     FUN= function(x) list(setNames(data.frame(data.ext[,x[1]]*data.ext[,x[2]]), 
                      paste(x, collapse="_")) )))

data.comb.sq <- do.call(cbind, sapply(input.vars,
                                       FUN = function(x) list(setNames(data.frame(data.ext[,x[1]]**2), 
                                                                      paste(x, "_sq", sep="")) )))


data.ext <- cbind(data.ext, data.comb.iter, data.comb.sq)


calc.model <- function(data, 
                       vars,
                       output,
                       base.model)
{
  
  
  data_sel <- data[,c(vars, output)]
  
  structure_base <- empty.graph(c(vars, output))
  modelstring(structure_base) <- base.model
  blacklist = data.frame(from = c(c("incidence", "death"),rep("incidence", length(vars)),rep("death", length(vars))),
                         to =  c(c("death", "incidence"),vars,vars))
  all_models <- list(base = list(struct = structure_base),
                     hc = list(struct = hc(data_sel, score = "bic-g", blacklist = blacklist))
                     #iamb_cor = list(struct = iamb(x = data_sel, test = "cor")),
                     #iamb_zf = list(struct = iamb(x = data_sel, test = "zf")),
                     #pc_cor = list(struct = pc.stable(x = data_sel, test = "cor")),
                     #pc_zf = list(struct = pc.stable(x = data_sel, test = "zf"))
  )
  for( name in names(all_models)){
    all_models[[name]][["fit"]] <- bn.fit(all_models[[name]][["struct"]], data = data_sel)
    all_models[[name]][["bic_score"]] <- bnlearn::score(all_models[[name]][["struct"]], data = data_sel, type = "bic-g")
    all_models[[name]][["arc_str"]] <- arc.strength(all_models[[name]][["struct"]], data = data_sel, criterion = "bic-g")
    all_models[[name]][["incidence_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$incidence$residuals)/var(data_sel$incidence))
    all_models[[name]][["death_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$death$residuals)/var(data_sel$death))
  }
  lm_models <- list(
    incidence = lm(incidence ~. ,data[c(vars,"incidence")]),
    death = lm(death ~. ,data[c(vars,"death")])
  )
  summary <- data.frame(model = c("lm", names(all_models)),
                        incidence_R2 = c(summary(lm_models$incidence)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["incidence_R2"]]})),
                        death_R2 = c(summary(lm_models$death)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["death_R2"]]}))
  )
  return(list(bn=all_models,
              lm=lm_models,
              summary=summary))
  
}


#mediation
models_H1 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand"),
                        c("incidence", "death"),
                        "[cases][vacc][HCA_supply][HCA_demand][incidence|cases:vacc:HCA_supply:HCA_demand][death|cases:vacc:HCA_supply:HCA_demand]")


#moderation
models_H2 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand",
                          "cases_HCA_supply", "cases_HCA_demand", "vacc_HCA_supply", "vacc_HCA_demand"),
                        c("incidence", "death"),
                        "[cases_HCA_supply][cases_HCA_demand][vacc_HCA_supply][vacc_HCA_demand][cases][vacc][HCA_supply][HCA_demand][incidence|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand][death|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand]")



models_H2$summary
models_H2$bn$hc
models_H2$bn$hc$arc_str
plot.network(models_H2$bn$hc)

models_H1$summary
models_H1$bn$hc
models_H1$bn$hc$arc_str

plot.network(models_H1$bn$hc)


#Step 1
fit.totaleffect = lm(incidence~cumulative_cases,data.all)
summary(fit.totaleffect) #0.27, R2 = 0.075
#Step 2
fit.mediator = lm(HCA_demand~cumulative_cases,data.all) 
summary(fit.mediator) #0.74, R2 = 0.55 
#Step 3

fit.dv = lm(incidence~HCA_demand+cumulative_cases,data.all)
summary(fit.dv) #0.36, 0.003 R2 = 0.13
library(tidyverse)
library(mediation)
#Step 4




results = mediate(fit.mediator, fit.dv, treat='cumulative_cases', mediator='HCA_demand', boot=T)



