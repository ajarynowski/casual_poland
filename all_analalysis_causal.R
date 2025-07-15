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
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}


set.seed(1020)

#https://www.medrxiv.org/content/10.1101/2021.12.08.21267167v2.full.pdf

data.all <- read.csv2("causal_data1.csv")

summary(data.all)

#[1] "X"                  "population_size"    "forest_density"     "PIS_support"        "income"             "postproduction_age"
#[7] "industry_revenue"   "empl_agriculture"   "poulation_density"  "Covid_deaths"       "size_COVID"         "betweenness_mob"   
#[13] "vacc"               "relig"              "HealthCareAcc"      "HealthAcc_phys"     "excess_mortality"   "inf_autumn"        
#[19] "tests"              "rt"                 "vacc_"              "county"             "deaths"             "deaths_norm"       
#[25] "HCA_demand"         "HCA_supply"         "cumulative_cases"   "vacc_fraq"          "deaths_"            "incidence"

names(data)

vars.sel.norm <- c("HCA_demand", "HCA_supply", 
              "cumulative_cases", "vacc_fraq", 
              "deaths_", "incidence")

vars.sel.raw <- c("HCA_demand", "HCA_supply", 
              "size_COVID", "vacc_fraq", 
              "population_size",
              "deaths", "inf_autumn")



#population
#0.032 56 942 biała podlaska 1.8
#0.170 296 401 białystok     1.8
#0.100 169 756 bielsko białą 1.69
#0.06 105 801 Koszalin       1.75
#0.14 244 104 Gdynia         1.74

hist(data.all$vacc_)

#data.all$population_size.fixed <- data.all$population_size*1780000



data.ext <- data.frame(cases = data.all$cumulative_cases,
                       vacc = data.all$vacc_fraq,
                       HCA_supply = data.all$HCA_supply,
                       HCA_demand = data.all$HCA_demand,
                       incidence = data.all$incidence,
                       death = ifelse(is.na(data.all$deaths_), 0, data.all$deaths_),
                       population_size = data.all$population_size,
                       poulation_density = 1.0*data.all$poulation_density,
                       income = 1.0*data.all$income)

input.vars <- c("cases", "vacc", "HCA_supply", "HCA_demand")

data.comb.iter <- do.call(cbind, combn(input_vars, 2, 
                     FUN= function(x) list(setNames(data.frame(data.ext[,x[1]]*data.ext[,x[2]]), 
                      paste(x, collapse="_")) )))

data.comb.sq <- do.call(cbind, sapply(input.vars,
                                       FUN = function(x) list(setNames(data.frame(data.ext[,x[1]]**2), 
                                                                      paste(x, "_sq", sep="")) )))


data.ext <- cbind(data.ext, data.comb.iter, data.comb.sq)



colnames(b)

data <- b
vars <- c("hits.fear","hits.mask","hits.quarantine")
output <- c("hits.covid-19","hits.covid-19")
base.model <- "[hits.fear][hits.mask][hits.quarantine][hits.covid-19|hits.quarantine:hits.mask:hits.fear][hits.covid-19:hits.covid-19|hits.quarantine:hits.mask:hits.fear]"

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
  }
  
  lm_models <- list(
    incidence = lm(incidence ~. ,data[c(vars,"hits.covid-19")])
  )
  
  summary <- data.frame(model = c("lm", names(all_models)),
                        incidence_R2 = c(summary(lm_models$incidence)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["incidence_R2"]]}))                        )
  
  return(list(bn=all_models,
              lm=lm_models,
              summary=summary))
}

models_H0 <- calc.model(data.ext,
                        c("cases", "vacc"),
                        c("incidence", "death"),
                        "[cases][vacc][incidence|cases:vacc][death|cases:vacc]")

models_H1 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand"),
                        c("incidence", "death"),
                        "[cases][vacc][HCA_supply][HCA_demand][incidence|cases:vacc:HCA_supply:HCA_demand][death|cases:vacc:HCA_supply:HCA_demand]")

models_H2 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand",
                          "cases_HCA_supply", "cases_HCA_demand", "vacc_HCA_supply", "vacc_HCA_demand"),
                        c("incidence", "death"),
                        "[cases_HCA_supply][cases_HCA_demand][vacc_HCA_supply][vacc_HCA_demand][cases][vacc][HCA_supply][HCA_demand][incidence|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand][death|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand]")



models_H3 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand",
                          "cases_HCA_supply", "cases_HCA_demand", "vacc_HCA_supply", "vacc_HCA_demand",
                          "population_size", "poulation_density", "cases_vacc"),
                        c("incidence", "death"),
                        "[cases_HCA_supply][cases_HCA_demand][vacc_HCA_supply][vacc_HCA_demand][cases][vacc][HCA_supply][HCA_demand][population_size][poulation_density][cases_vacc][incidence|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand:population_size:poulation_density:cases_vacc][death|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand:population_size:poulation_density:cases_vacc]")


#names(data.comb.iter)

models_H4 <- calc.model(data.ext,
                        c("cases", "vacc", "HCA_supply", "HCA_demand",
                          "cases_HCA_supply", "cases_HCA_demand", "vacc_HCA_supply", "vacc_HCA_demand",
                          
                          "population_size", "poulation_density",
                          "cases_sq", "vacc_sq", "HCA_supply_sq", "HCA_demand_sq", 
                          "cases_vacc", "HCA_supply_HCA_demand"
                          ),
                        c("incidence", "death"),
                        "[cases_HCA_supply][cases_HCA_demand][vacc_HCA_supply][vacc_HCA_demand][cases][vacc][HCA_supply][HCA_demand][population_size][poulation_density][cases_sq][vacc_sq][HCA_supply_sq][HCA_demand_sq][cases_vacc][HCA_supply_HCA_demand][incidence|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand:population_size:poulation_density:cases_sq:vacc_sq:HCA_supply_sq:HCA_demand_sq:cases_vacc:HCA_supply_HCA_demand][death|cases:vacc:HCA_supply:HCA_demand:cases_HCA_supply:cases_HCA_demand:vacc_HCA_supply:vacc_HCA_demand:population_size:poulation_density:cases_sq:vacc_sq:HCA_supply_sq:HCA_demand_sq:cases_vacc:HCA_supply_HCA_demand]")



models_H0$summary
models_H1$summary
models_H2$summary
models_H3$summary
models_H4$summary
models_H1$bn$hc
models_H2$bn$hc
models_H3$bn$hc
models_H4$bn$hc


plot.network(models_H0$bn$hc$struct)
plot.network(models_H1$bn$hc$struct)
plot.network(models_H2$bn$hc$struct)
plot.network(models_H3$bn$hc$struct)
plot.network(models_H4$bn$hc$struct)



models_H0$summary
models_H0$bn$hc$arc_str
models_H1$summary
models_H1$bn$hc$arc_str
models_H2$summary
models_H2$bn$hc$arc_str
models_H3$summary
models_H3$bn$hc$arc_str
models_H4$summary
models_H4$bn$hc$arc_str

models_H2$bn$hc
models_H3$bn$hc
models_H4$bn$hc


#model without normalization

load("infekcje_2022_caly.RData")
inf_2022_all$deaths = ifelse(is.na(inf_2022_all$deaths), 0, inf_2022_all$deaths)

data <- inf_2022_all
vars <- c("size_COVID", "vacc")
output <- c("inf_autumn", "deaths")
base.model <- "[size_COVID][vacc][inf_autumn|size_COVID:vacc][death|size_COVID:vacc]"

calc.model <- function(data, 
                       vars,
                       output,
                       base.model)
{
  

data_sel <- data[,c(vars, output)]

structure_base <- empty.graph(c(vars, output))
modelstring(structure_base) <- base.model
blacklist = data.frame(from = c(c("inf_autumn", "deaths"),rep("inf_autumn", length(vars)),rep("deaths", length(vars))),
                       to =  c(c("deaths", "inf_autumn"),vars,vars))
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
  all_models[[name]][["incidence_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$inf_autumn$residuals)/var(data_sel$inf_autumn))
  all_models[[name]][["death_R2"]] <- 1 - (var(all_models[[name]][["fit"]]$deaths$residuals)/var(data_sel$deaths))
}
lm_models <- list(
  incidence = lm(inf_autumn ~. ,data[c(vars,"inf_autumn")]),
  death = lm(deaths ~. ,data[c(vars,"deaths")])
)
summary <- data.frame(model = c("lm", names(all_models)),
                      incidence_R2 = c(summary(lm_models$incidence)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["incidence_R2"]]})),
                      death_R2 = c(summary(lm_models$death)$r.squared, sapply(names(all_models), FUN = function(name){all_models[[name]][["death_R2"]]}))
)
return(list(bn=all_models,
            lm=lm_models,
            summary=summary))

}
models_H1 <- calc.model(inf_2022_all,
                        c("size_COVID", "vacc", "HCA_supply", "HCA_demand"),
                        c("inf_autumn", "deaths"),
                        "[size_COVID][vacc][HCA_supply][HCA_demand][inf_autumn|size_COVID:vacc:HCA_supply:HCA_demand][deaths|size_COVID:vacc:HCA_supply:HCA_demand]")

models_H1$summary
plot.network(models_H1$bn$hc$struct)


models_H0 <- calc.model(inf_2022_all,
                        c("size_COVID", "vacc"),
                        c("inf_autumn", "deaths"),
                        "[size_COVID][vacc][inf_autumn|size_COVID:vacc][deaths|size_COVID:vacc]")


models_H0$summary
plot.network(models_H0$bn$hc$struct)



models_H1 <- calc.model(inf_2022_all,
                        c("size_COVID", "vacc", "HCA_supply", "HCA_demand"),
                        c("inf_autumn", "deaths"),
                        "[size_COVID][vacc][HCA_supply][HCA_demand][inf_autumn|size_COVID:vacc:HCA_supply:HCA_demand][deaths|size_COVID:vacc:HCA_supply:HCA_demand]")
models_H1$summary
plot.network(models_H1$bn$hc$struct)



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

ggplot (iot3, aes(date, as.numeric(hits))) + 
  geom_line(aes(color=keyword))+
  theme_bw()

results = mediate(fit.mediator, fit.dv, treat='Sepal.Length', mediator='mediator', boot=T)

results = mediate(fit.mediator, fit.dv, treat='cumulative_cases', mediator='HCA_demand', boot=T)



ccfvalues = ccf(b$hits.fear,b$`hits.covid-19`, na.action = na.pass, lag.max=7)
