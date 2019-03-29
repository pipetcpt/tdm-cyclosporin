# setup ----

calculate_clearance <- function(WEIGHT, AGE, TBIL, POD_week){
  TVCL <- 28.5 - 1.24*POD_week - 0.252*(TBIL-10) + 0.188*(WEIGHT-60) - 0.191*(AGE - 40)
}

TVV1 <- 133.00
TVKA <-   1.28
    
default_dose_example_csv <- '"Date","Dosing_Time","Route","Dose"
"17.05.03","10:30","PO","750"
"17.05.04","10:30","PO","750"
"17.05.04","22:30","PO","500"'

input <- list( # These values should be defined in ui.R.
  sex = 'male',
  obsc = 1000, # ng/mL 
  obsDate = "2017-05-05", 
  obsTime = strptime("10:30", "%R"),
  weight = 60,
  Observations = '1',
  total_bilirubin = 1,
  scr = 0.9,
  post_op_date = 3,
  newdose = 1000,
  newtau = 48,
  newinf = 1,
  ll = 15,
  ul = 40,
  age = 20
)

library(deSolve)
library(plyr)
library(grid)
library(compiler)
library(shinyTime)
library(lubridate)
library(TeachingDemos)
library(rmarkdown)
library(knitr)
library(DT)
library(rsconnect)
library(tidyverse)

# dose in mg kg당 6-8mg : 8*60

## ltv2mat copy right:: Prof. Bae 

cmat=function(vec){
  LENGTH=length(vec)
  DIM=as.integer(round((sqrt(8*LENGTH+1)-1)/2,0))
  if(DIM*(DIM+1)/2!=LENGTH) return(NULL)
  mat=matrix(nrow=DIM, ncol=DIM)
  mat[upper.tri(mat, diag=TRUE)]=vec
  mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
  return(mat)
}

# main ----

shiny::shinyServer(function(input, output) {
  
  # chapter 1. downloadData ----
  
  output$downloadData <- downloadHandler(
    filename <- function() { paste("dose_example", '.csv', sep='') },
    content <- function(file) {
      write.csv(datasetInput(), file, row.names = F)
    }
  )
  
  # Chapter 2. dosing_history_contents ----
  output$dosing_history_contents <- renderTable({
    if (is.null(input$file1)) return(read_csv(default_dose_example_csv,
                                              col_types = 'cccd'))
    return(read_csv(input$file1$datapath,
                    col_types = 'cccd'))
  })
  
  # Chapter 3. creatinine_clearance ----
  # crcl: https://www.mdcalc.com/creatinine-clearance-cockcroft-gault-equation
  output$typical_values <- renderTable({
    # Typical Values of Population PK parameters
    TVCL <- calculate_clearance(input$weight,input$age,input$total_bilirubin,input$post_op_date/7) %>% 
      round(digits = 2)
    return(tibble(`CL/F (L/h)` = TVCL, `V/F (L)` = TVV1, `Ka (/hr)` = TVKA))
  })
  
  # Chapter 4. output_table1_time_predicted_concentration ----
  output$output_table1_time_predicted_concentration <- renderTable({
    sim_data_output <- sim.data()
    prtx_predicted_concentration <- sim_data_output$table1
    return(prtx_predicted_concentration)
  })
  
  # Chapter 5. outputtable2 ----
  output$outputtable2 <- renderTable({
    sim_data_output <- sim.data()
    prtx_predicted_concentration <- sim_data_output$table2
    return(prtx_predicted_concentration)
  })
  
  # Chapter 6. outputtable3 ----
  output$outputtable3 <- renderTable({
    prt1=sim.data()
    prt2=prt1[complete.cases(prt1),]
    if(input$Observations=='1'){
    }
    if (input$Observations=='2')
    {
      prtx=prt2[2,c("CL","V1","V2")]
      return(prtx)
    }
  })
  
  # prelude 1. datasetInput ----
  
  datasetInput <- reactive({
    coln=c("Date","Inf_st_Time","Inf_ed_Time","Dose" )
    dat1=c("17.05.03", "10:30", "11:30","500")
    dat2=c("17.05.03", "22:30", "23:30","500")
    dat3=c("17.05.04", "10:30", "11:30","1000")
    dat4=c("17.05.04", "22:30", "23:30","1000")
    
    suppl <- data.frame (rbind(dat1,dat2,dat3,dat4))
    colnames(suppl)<-coln
    rownames(suppl)<-NULL
    suppl
  })
  
  # prelude 2. dose.data ----
  
  dose.data <- reactive({
    inFile <- input$file1 
    if (is.null(inFile)) return(NULL) 
    a=read.csv(inFile$datapath, header=T, stringsAsFactors = T) 
    b=a[complete.cases(a), ]  
    b$paste=paste(b$Date,b$Time) 
    b
    #Time calculation code is copyrighted
  })
  
  sim.data <- reactive({
    
    # prelude 3. sim.data ----
    
    # input: Observation
    obs1conc <- input$obsc
    obs1time <- input$obst
    obs1dat  <- input$obsd
    
    # input: Demog
    
    # Typical Values
    
    TVCL <- calculate_clearance(input$weight,input$age,input$total_bilirubin,input$post_op_date/7)
    
    # Eta (Omega)
    ETA1SD <- 0.04     # Vancomycin Clearance eta 
    ETA2SD <- 0.04     # Vancomycin Volume eta
    omega <- cmat(c(ETA1SD ,0,ETA2SD)) %>%  # Prof. Bae's function
      print() #c(Clearance eta,0,volume eta)
    omega.inv <- solve(omega) %>% 
      print()
    
    # Eps (Sigma)
    EPS1SD <- 40      # Vancomycin Additive residual error   
    EPS2SD <- 0.1     # Vancomycin Proportional residual error
    EPS2SDsq=(EPS2SD)^2 # Vancomycin square
    
    input_file_text <- ifelse(is.null(input$file1), 
                              yes = default_dose_example_csv, 
                              no = input$file1$datapath)
    
    rawdata <- read_csv(input_file_text, col_types = 'cccd') %>% 
      as.data.frame() %>% 
      print()
    
    dosing_data <- rawdata %>% 
      mutate(date_time = sprintf('%sT%s', Date, Dosing_Time)) %>% 
      mutate(date_time_lubridate = ymd_hm(date_time)) %>% 
      mutate(Concentration = NA) %>% 
      select(date_time_lubridate, Dose, Concentration) %>% 
      as_tibble() %>% 
      print()
    
    observation_data <- tibble(date_time_lubridate = ymd_hms(sprintf('%sT%s', 
                                                                     input$obsDate, 
                                                                     substr(input$obsTime, 12, 20))),
                               Dose = NA,
                               Concentration = input$obsc) %>% 
      print()
    
    Dose_Concentration <- bind_rows(dosing_data, observation_data) %>% 
      mutate(first_dose = first(date_time_lubridate)) %>% 
      mutate(actual_time = difftime(date_time_lubridate, first_dose, unit = 'hour') %>% 
               as.numeric() %>% 
               round(2)) %>% 
      print()
    
    DOSEdata<- Dose_Concentration %>% 
      filter(!is.na(Dose)) %>% 
      mutate(var = 1, method = 'add') %>% 
      select(var, time = actual_time, value = Dose, method) %>% 
      as.data.frame() %>% 
      print()
    
    
    
    # model ----
    
    model <- function(Time, A, eta){
      Ka <- TVKA
      Cl <- TVCL * exp(eta[1])
      Vd <- TVV1 * exp(eta[2])
      
      Ke <-  Cl/Vd  # Elimination
      
      dA <- vector(length = 2)
      dA[1] <- -Ka * A[1]               # Depot compartment - intestine 
      dA[2] <-  Ka * A[1] - Ke * A[2]   # Central compartment - plasma
      return(list(dA))
    }
    
    mod.cmp <-  compiler::cmpfun(model)
    
    initial_amount <- c(A1 = 0, A2 = 0) 
    
    Concentration_only <- Dose_Concentration %>% 
      filter(!is.na(Concentration)) %>% 
      print()
    Dose_only <- Dose_Concentration %>% 
      filter(!is.na(Dose)) %>% 
      print()
    
    
    Observeddate <- Concentration_only$date_time_lubridate
    pointtime <- Concentration_only$actual_time
    TIME <- seq(from = 0, to = pointtime, by = 0.1)
    TIMElast <- max(TIME)
    
    if (input$Observations=='1') {
      y <- obs1conc / 1000
      #' @example mapb2(c(0.1, 0.1))
      mapb2 <- function(eta){ # eta is a list of 2. eta <- c(0.04321,0.04321)
        etamat=matrix(unlist(eta))
        out <- deSolve::lsoda(y = initial_amount, 
                              times = TIME, 
                              func = model, 
                              parms = eta, 
                              events=list(data=DOSEdata)) %>% 
          as_tibble() %>% 
          mutate(DV = A2/TVV1)
        
        eta <- c(eta[1],eta[2])
        eta_m <- unlist(matrix(eta,nrow = 2))
        sig2 <- EPS2SDsq
        sig2j <- subset(out[,4],out[,1]==pointtime)^2*sig2
        sqwres <- log(sig2j) + (1/sig2j) * (y[1]-subset(out[,4],out[,1]==pointtime))^2
        nOn <- diag(t(eta_m) %*% omega.inv %*% eta_m)
        return(sum(sqwres)+ nOn)
      } 
    }
    
    mapb2.cmp <- cmpfun(mapb2)
    ini <- c(0.04,0.04)
    mapb2.cmp(ini)
    
    # Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound. The initial value must satisfy the constraints. This uses a limited-memory modification of the BFGS quasi-Newton method. If non-trivial bounds are supplied, this method will be selected, with a warning.
    
    shiny::withProgress(
      message = 'Minimization in progress', 
      min = 0, 
      max = 100,
      value = 99, 
      {
        FIT <- stats::optim(par = ini, # CL, V => fitting
                            fn = mapb2.cmp, # A function to be minimized (or maximized)
                            method="L-BFGS-B", 
                            control = list(trace=TRUE,REPORT=TRUE))
        print(FIT$par)
        
      }
    )
    
    outs <- lsoda(y = initial_amount, 
                  times = TIME, 
                  func = mod.cmp,  
                  parms = FIT$par,
                  events = list(data = DOSEdata)) %>% 
      as.data.frame() %>% 
      mutate(DV = A2/TVV1) %>% # mg / L == ug/mL == ng / uL
      print()
    
    if (input$Observations=='1'){
      outs$pointtime=pointtime  
      outs$Observeddate=Observeddate
      outs$predictedConc=subset(outs[,4],outs[,1]==pointtime)
      outs$observedConc=obs1conc
    }
    
    outs$TIME <- TIME
    outs$CL <- TVCL*(exp(FIT$par[1]))    #predicted CL
    outs$V1 <- TVV1*(exp(FIT$par[2]))    #predicted V2
    
    outs2=merge(x=outs,y=DOSEdata, by="time",all.x=TRUE)
    head(outs2);tail(outs2)
    outs2 %>% as_tibble()
    
    table1_raw <- tibble(time = pointtime,
                         obs_conc = outs %>% 
                           filter(time == pointtime) %>% 
                           .$observedConc,
                         pred_conc = outs %>%
                           filter(time == pointtime) %>% 
                           .$predictedConc * 1000)
    
    cyclosporine_pk_plot <- outs2 %>% 
      as_tibble() %>% 
      ggplot(aes(time, DV*1000)) + 
      geom_point(data = table1_raw, aes(x=pointtime, y=pred_conc, color = 'Predicted')) + 
      geom_point(data = table1_raw, aes(x=pointtime, y=obs_conc, color = 'Observed')) + 
      geom_line(alpha = 0.5) +
      geom_hline(yintercept = c(50, 200), color = 'red') +
      labs(x = 'Time (hour)', y = 'Cyclosporine concentration (ng/mL)',
           title = 'Oral administration of cyclosporine - surgical use',
           subtitle = 'Therapeutic Drug Monitoring',
           color = 'Cyclosporine (ng/mL)') +
      theme_bw()
    cyclosporine_pk_plot
    
    sim_data_output <- list(
      FIT = FIT,
      table1 = table1_raw %>% select(`time (h)` = 1,
                                     `observed conc. (mg/L)` = 2,
                                     `predicted conc. (mg/L)` = 3),
      table2 = tibble(`CL/F (L/h)` = TVCL*(exp(FIT$par[1])),
                      `V/F (L)` = TVV1*(exp(FIT$par[2]))),
      plot1 = cyclosporine_pk_plot,
      DOSEdata = DOSEdata,
      initial_amount = initial_amount,
      mod.cmp = mod.cmp)
    
    return(sim_data_output)
  })
  
  output$plotCONC <- renderPlot({
    sim_data_output <- sim.data()
    return(sim_data_output$plot1)
  })
  
  output$plotCONC2 <- renderPlot({
    sim_data_output <- sim.data()
    
    DOSEdata <- sim_data_output$DOSEdata
    initial_amount <- sim_data_output$initial_amount
    mod.cmp <- sim_data_output$mod.cmp
    FIT <- sim_data_output$FIT
    
    last_dose_actual_time <- max(DOSEdata$time)
    future_tau  <- input$newtau
    future_dose <- input$newdose
    
    DOSEdata_future <- bind_rows(DOSEdata,  
                                 tibble(seq = 1:10,
                                        actual_time = last_dose_actual_time + seq * future_tau,
                                        Dose = future_dose) %>% 
                                   mutate(var = 1, method = 'add') %>% 
                                   select(var, time = actual_time, value = Dose, method)) %>% 
      print()
    
    TIME_future <- seq(from = 0, to = max(DOSEdata_future$time) , by = 0.1)
    outs_future <- lsoda(y = initial_amount, 
                         times = TIME_future, 
                         func = mod.cmp,  
                         parms = FIT$par,
                         events = list(data = DOSEdata_future)) %>% 
      as.data.frame() %>% 
      mutate(DV = A2/TVV1) %>% # mg / L == ug/mL == ng / uL
      as_tibble() %>% 
      print()
    
    cyclosporine_dose_adjustment_plot <- outs_future %>% 
      as_tibble() %>% 
      ggplot(aes(time, DV*1000)) + 
      #geom_point(data = table1_raw, aes(x=pointtime, y=pred_conc, color = 'Predicted')) + 
      #geom_point(data = table1_raw, aes(x=pointtime, y=obs_conc, color = 'Observed')) + 
      geom_line(alpha = 0.5) +
      geom_hline(yintercept = input$ul, color = 'red') +
      geom_hline(yintercept = input$ll, color = 'blue') +
      labs(x = 'Time (hour)', y = 'Cyclosporine concentration (ng/mL)',
           title = 'Oral administration of cyclosporine - surgical use',
           subtitle = 'Therapeutic Drug Monitoring',
           color = 'Cyclosporine (ng/mL)') +
      scale_x_continuous(limits = c(max(DOSEdata$time), max(DOSEdata_future$time)))+
      theme_bw()
    cyclosporine_dose_adjustment_plot
    return(cyclosporine_dose_adjustment_plot)
  })
  
  
  # end ----  
})
