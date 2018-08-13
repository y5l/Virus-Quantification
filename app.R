#source("https://bioconductor.org/biocLite.R")
#biocLite("flowCore")

#bringing in several packages that will give useful functions for the app
library(ggplot2)
library(flowCore)
library(stringr)
library(dplyr)
library(reshape2)
library(shiny)
library(shinyFiles)


#note: in the comments, block number and plate location mean the same thing

####################################################################### User Interface portion of app ############################################################################

                                                                                                   # Defines UI for application (must be stored into variable ui), fluidPage makes it so the page consists of rows and columns and automatically scales to the browser size
ui <- fluidPage(
   
                                                                                                   # titlePanel takes a string and sets it as the title of the page
   titlePanel("Titer Estimation from Flow Cytometer Side Scatter"),
   
                                                                                                   # sidebarPanel takes some arguments and adds them to a side bar panel
   sidebarPanel(
     
                                                                                                   # fileInput creates a box which allows adding a file, inputID gives the inputted data an ID x and effectively stores it as a column of a variable input, so it can be accessed using input$x, label is just the title of the file input box, and multiple = TRUE allows the selection of more than 1 file
     fileInput(inputId = "files",label = "Choose file",multiple = TRUE),                                                                                          
     
                                                                                                   # numericInput creates a box which allows the input of a number, inputId and label have the same effect as the case with fileInput, min and max are the minimum and maximum values of the allowed inputs, and value is the default value of the input
     numericInput(inputId = "FSCgatelower",label = "Lower FSC Gate (Should be numeric and less than upper FSC gate)", min = 0, max = 50e6, value = 3.125e6),
     numericInput(inputId = "FSCgateupper",label = "Upper FSC Gate (Should be numeric and greater than lower FSC gate)", min = 0, max = 50e6, value = 14e6),
     numericInput(inputId = "ssclowlimit", label = "Lower Limit for Calculating Percent SSC", min = 0, max = 1e7, value = 1e6)
   ),
   
                                                                                                   # mainPanel takes in several inputs and adds them to the main panel
   mainPanel(
     
                                                                                                   # tabsetPanel creates a panel that contains a set of tabs, type can either be tab or pills, pills looks likes tab but is more simple
     tabsetPanel(type = "tab",
                 
                                                                                                   # tabPanel creates a new tab in the tabsetPanel with its title being title, followed by output elements to include in the tab
                 tabPanel(title = "Estimated Titers", tableOutput(outputId = 'dataOutput')),  
                 tabPanel(title = "Calibration Curve", plotOutput(outputId = 'curveOutput')),  
                 tabPanel(title = "Percent Side Scatter ~ Dilution", plotOutput(outputId = 'persscdilOutput')), 
                 tabPanel(title = "SSC ~ FSC", plotOutput(outputId = 'fscsscOutput')),  
                 tabPanel(title = "Negative StDev, LOD, LOQ", tableOutput(outputId = 'negsdlodloqOutput')),
                 tabPanel(title = "FSC and SSC histograms", plotOutput(outputId = 'fscsschistOutput'),numericInput(inputId = "histbins", label= "Number of bins for the FSC and SSC histograms", min=1,max=1000,value=30)),
                 tabPanel(title = "SSC and Percent SSC for Blank Blocks", plotOutput(outputId = 'sscblankplotOutput'),tableOutput(outputId = "sscblanktableOutput"))
     )
  
   )
)

######################################################################## Server portion of app ####################################################################################

                                                                                                   # Defines what is being done at the server for an application (must be function(input,output) stored into variable server), 
server <- function(input, output) {

############################################################################### Inputs ############################################################################################   
  
                                                                                                   # Making a "reactive" variable (which will contain a table, coded for later) that changes as the files that are inputted change
  dataInput <- reactive({
    
                                                                                                   # input$files was defined when we set the inputId of the fileinput to be "files", this then stores the type of file into a variable Type
      Type <- input$files$type
                                                                                                   # initiating 2 more variables files and name as empty vectors
      files <- vector()
      name <- vector()
                                                                                                   # going through each entry of Type, if it isn't a csv file, its data path is added to the files variable and the name of the file is added onto names
      for (i in 1:length(Type)) {
        if (Type[i] !="text/csv"){
          files <- c(files,input$files$datapath[i])
          name <- c(name,input$files$name[i])
        }
      }

                                                                                                   # making a data frame with both the name and files vectors (we will need this later)
     name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
     
                                                                                                   # going through each file again, but this time, if it IS a csv file, store the data path into csvtext
     for (i in 1:length(Type)){
        if (Type[i] =="text/csv"){
          csvtext<- input$files$datapath[i]
        }
     }
     
                                                                                                   # stores the actual data from csvtext (which is a data path), into METADat
     METADat <- read.csv(csvtext)
     
                                                                                                   # initiates variable O_Meta_Dat (Overall Meta Data), which will be used later
     O_Met_Dat <- 0
    
                                                                                                   # making a function overall that takes in the file path
     createoverall <- function(filepath) {
       
                                                                                                   # this block stores the SSC and FSC into a data frame impDAT and then gives them the respective column names
       flow <- read.FCS(file.path(filepath))

       impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
       colnames(impDAT) <- c('SSC','FSC')
       
                                                                                                   # getting the actual file name from the file path, by using the name.and.files dataframe created earlier
       filename <- name.and.files$name[name.and.files$files==filepath]

                                                                                                   # modifying the file name to get the well, row, and column
       filename.split<- gsub('.fcs','',filename)
       well<-filename.split
       row<-gsub("[[:digit:]]", "", filename.split)
       column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))
       
                                                                                                   # makes sure the FSC column of impDAT is numeric
       impDAT$FSC<-as.numeric(as.character(impDAT$FSC))
       
       
                                                                                                   # filters the data from impDAT with an upper and lower FSC gate, storing the result to modDAT
       modDAT<-impDAT[impDAT$FSC>=input$FSCgatelower & impDAT$FSC<=input$FSCgateupper,]
       
                                                                                                   # gets the number of cells from the FSC column of impDAT
       Number_of_Cells<-length(impDAT$FSC)
       
                                                                                                   # stores into Percent.FSC the percent of entries in the FSC column of impDAT that are present in that of modDAT
       Percent.FSC<-length(modDAT$FSC)/length(impDAT$FSC)*100
       
                                                                                                   # stores into Percent.SSC the percent of entries in the modDAT$FSC column that are in the modDAT$SSC column and are greater than the SSC lower limit
       Percent.SSC<-length(which(as.numeric(as.character(modDAT$SSC))>input$ssclowlimit))/length(modDAT$FSC)*100

                                                                                                   # there might be a more succinct way of programming this section
                                                                                                   # these if statements check the position of the well, and stores data associated with that well into several variables
       if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==1]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==1]
         Type<-METADat$Type[METADat$Plate_Location==1]
         Plate_Location <- 1
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==1]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==1]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==1]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==1]}
       }
       
       if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==2]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==2]
         Type<-METADat$Type[METADat$Plate_Location==2]
         Plate_Location <- 2
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==2]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==2]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==2]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==2]}
       }
       
       if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==3]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==3]
         Type<-METADat$Type[METADat$Plate_Location==3]
         Plate_Location <- 3
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==3]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==3]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==3]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==3]}
       }
       
       if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==4]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==4]
         Type<-METADat$Type[METADat$Plate_Location==4]
         Plate_Location <- 4
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==4]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==4]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==4]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==4]}
       }
       
       if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==5]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==5]
         Type<-METADat$Type[METADat$Plate_Location==5]
         Plate_Location <- 5
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==5]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==5]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==5]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==5]}
       }
       
       if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
         ID<- METADat$ID[METADat$Plate_Location==6]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==6]
         Type<-METADat$Type[METADat$Plate_Location==6]
         Plate_Location <- 6
         if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==6]}
         if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==6]}
         if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==6]}
         if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==6]}
       }
       
       if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==7]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==7]
         Type<-METADat$Type[METADat$Plate_Location==7]
         Plate_Location <- 7
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==7]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==7]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==7]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==7]}
       }
       
       if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==8]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==8]
         Type<-METADat$Type[METADat$Plate_Location==8]
         Plate_Location <- 8
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==8]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==8]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==8]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==8]}
       }
       
       if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==9]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==9]
         Type<-METADat$Type[METADat$Plate_Location==9]
         Plate_Location <- 9
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==9]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==9]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==9]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==9]}
       }
       
       if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==10]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==10]
         Type<-METADat$Type[METADat$Plate_Location==10]
         Plate_Location <- 10
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==10]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==10]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==10]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==10]}
       }
       
       if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==11]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==11]
         Type<-METADat$Type[METADat$Plate_Location==11]
         Plate_Location <- 11
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==11]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==11]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==11]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==11]}
       }
       
       if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
         ID<- METADat$ID[METADat$Plate_Location==12]
         PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==12]
         Type<-METADat$Type[METADat$Plate_Location==12]
         Plate_Location <- 12
         if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==12]}
         if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==12]}
         if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==12]}
         if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==12]}
       }
       
                                                                                                   # stores all the variables into a data.frame overall, for the current file path
       overall <- data.frame(filename = filename, Type= Type, well = well, row = row,
                             column = column, Dilution = Dilution, ID = ID, PA_Titer = PA_Titer,
                             Percent.SSC = Percent.SSC, Percent.FSC= Percent.FSC, Number_of_Cells=Number_of_Cells, Plate_Location = Plate_Location)
                                                                  # returns the overall dataframe
       return(overall)
     }
     
                                                                                                   # applies the createoverall function to each item in the files list (which is a list of data paths)
     z <- lapply(files,createoverall)
     
                                                                                                   # z is a vertical list of horizontal lists, so it needs to be reformatted into a dataframe, which the line below does (there is likely a better way to program this)
     O_Met_Dat <- rbind(z[[1]],z[[2]],z[[3]],z[[4]],z[[5]],z[[6]],z[[7]],z[[8]],z[[9]],z[[10]],z[[11]],z[[12]],z[[13]],z[[14]],z[[15]],z[[16]],z[[17]],z[[18]],z[[19]],z[[20]],z[[21]],z[[22]],z[[23]],z[[24]],z[[25]],z[[26]],z[[27]],z[[28]],z[[29]],z[[30]],z[[31]],z[[32]],z[[33]],z[[34]],z[[35]],z[[36]],z[[37]],z[[38]],z[[39]],z[[40]],z[[41]],z[[42]],z[[43]],z[[44]],z[[45]],z[[46]],z[[47]],z[[48]],z[[49]],z[[50]],z[[51]],z[[52]],z[[53]],z[[54]],z[[55]],z[[56]],z[[57]],z[[58]],z[[59]],z[[60]],z[[61]],z[[62]],z[[63]],z[[64]],z[[65]],z[[66]],z[[67]],z[[68]],z[[69]],z[[70]],z[[71]],z[[72]],z[[73]],z[[74]],z[[75]],z[[76]],z[[77]],z[[78]],z[[79]],z[[80]],z[[81]],z[[82]],z[[83]],z[[84]],z[[85]],z[[86]],z[[87]],z[[88]],z[[89]],z[[90]],z[[91]],z[[92]],z[[93]],z[[94]],z[[95]],z[[96]])

                                                                                                   # makes sure the PA_Titer and Dilution columns are numeric, and also adds Act.Tit (Actual titer) and LN_Titer (natural log Titer) columns
     O_Met_Dat$PA_Titer<-as.numeric(as.character(O_Met_Dat$PA_Titer))
     O_Met_Dat$Dilution<-as.numeric(as.character(O_Met_Dat$Dilution))
     O_Met_Dat$Act.Tit<-O_Met_Dat$PA_Titer/O_Met_Dat$Dilution
     O_Met_Dat$LN_Titer<-log(O_Met_Dat$Act.Tit)
     
                                                                                                   # gets the mean and standard deviation of the SSC values whose type is blank, and stores that to Negative and NegStDev, respectively
     Negative<-mean(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
     NegStDev<-sd(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
     
                                                                                                   # limit of detection is 3 times the NegStDev, limit of quantification is 10 times
     LOD<-3*NegStDev
     LOQ<-10*NegStDev
 
                                                                                                   # an Adj.Per (Adjusted percentage) column is added by subtracting Negative from the Percent.SSC values
     O_Met_Dat$Adj.Per<-O_Met_Dat$Percent.SSC-Negative
     
 
                                                                                                   # a calibration curve is created from the LN_Titer and Adj.Per columns, from both only which the values with Type == Standard and the corresponding Adj.Per value is greater than the LOQ are considered
     Calib<-lm(O_Met_Dat$LN_Titer[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ]~O_Met_Dat$Adj.Per[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ])    
     
                                                                                                   # titer estimations are calculated using the calibration curve and values from the Adj.Per column, and stored into the FC_Titer column in rows for which the Adj.Per value is greater than the LOQ
     O_Met_Dat$FC_Titer[O_Met_Dat$Adj.Per>LOQ] <-exp(Calib$coefficients[1] + O_Met_Dat$Adj.Per[O_Met_Dat$Adj.Per>LOQ] * Calib$coefficients[2]) * O_Met_Dat$Dilution[O_Met_Dat$Adj.Per>LOQ]
                                                                                                   # calculates the mean titer for each block
     titers <- vector()
     totalstatus <- vector()
     
                                                                                                   #going through each block number
     for (i in 1:12) {
                                
                                                                                                   #initiates 2 vectors to store the titer estimation of each block, and the status of the block (in terms of being above or below the LOD and LOQ)
       location <- vector()
       blockstatus <- vector()
       for (j in 1:96) {
                                                                                                   # status = 0 assumes that the adjusted percentage is above the LOQ, = 1 assumes b/w LOQ and LOD, = 2 assumes below LOD
         status <- 0
         if (O_Met_Dat$Plate_Location[j] == i) {
           location <- c(location, O_Met_Dat$FC_Titer[j])
           if (O_Met_Dat$Adj.Per[j] < LOQ) {
             status <- 1
           }
           if (O_Met_Dat$Adj.Per[j] < LOD) {
             status <- 2
           }
           
                                                                                                   # adds the status values to a vector so the statuses of the block can be evaluated later
           blockstatus <- c(blockstatus,status)
         }
       }
       
                                                                                                   #assigns a message to the block, into totalstatus, based on the highest status value in the status block
       if (max(blockstatus)==0) {
         totalstatus[i] <- "All Adjusted Percentages are above LOQ"
       }
       if (max(blockstatus)==1) {
         totalstatus[i] <- "Some Adjusted Percentages are below LOQ but above LOD"
       }
       if (max(blockstatus)==2) {
         totalstatus[i] <- "Some Adjusted Percentages are below LOD"
       }
        
                                                                                                   # takes the mean of the titers from each well and stores it into titers
       titers <- c(titers,mean(location))
     }
                                                                                              # formats the table into scientific notation and sets its dimension to be 6x2
     titers <- formatC(titers, format = "e", digits = 2)

                                                                                       # formatting the data frames and adding a plate location column so the table is ready to be displayed
     pl <- as.data.frame(c(1:12))
     titers <- cbind(pl,titers,totalstatus)
     colnames(titers) <- c("Plate Location","Titer","Status")
     
                                                                                                   # calling up titers so the data within gets displayed
     titers
   })
  
  curveInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    for (i in 1:length(Type)){
      if (Type[i] =="text/csv"){
        csvtext<- input$files$datapath[i]
      }
    }
    METADat <- read.csv(csvtext)
    O_Met_Dat <- 0
    
    createoverall <- function(filepath) {
      flow <- read.FCS(file.path(filepath))
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filepath]
      
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      row<-gsub("[[:digit:]]", "", filename.split)
      column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))

      impDAT$FSC<-as.numeric(as.character(impDAT$FSC))
      
      modDAT<-impDAT[impDAT$FSC>=input$FSCgatelower & impDAT$FSC<=input$FSCgateupper,]
      Number_of_Cells<-length(impDAT$FSC)
      Percent.FSC<-length(modDAT$FSC)/length(impDAT$FSC)*100
      Percent.SSC<-length(which(as.numeric(as.character(modDAT$SSC))>input$ssclowlimit))/length(modDAT$FSC)*100

      if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==1]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==1]
        Type<-METADat$Type[METADat$Plate_Location==1]
        Plate_Location <- 1
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==1]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==1]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==1]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==1]}
      }
      
      if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==2]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==2]
        Type<-METADat$Type[METADat$Plate_Location==2]
        Plate_Location <- 2
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==2]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==2]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==2]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==2]}
      }
      
      if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==3]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==3]
        Type<-METADat$Type[METADat$Plate_Location==3]
        Plate_Location <- 3
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==3]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==3]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==3]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==3]}
      }
      
      if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==4]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==4]
        Type<-METADat$Type[METADat$Plate_Location==4]
        Plate_Location <- 4
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==4]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==4]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==4]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==4]}
      }
      
      if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==5]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==5]
        Type<-METADat$Type[METADat$Plate_Location==5]
        Plate_Location <- 5
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==5]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==5]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==5]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==5]}
      }
      
      if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==6]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==6]
        Type<-METADat$Type[METADat$Plate_Location==6]
        Plate_Location <- 6
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==6]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==6]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==6]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==6]}
      }
      
      if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==7]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==7]
        Type<-METADat$Type[METADat$Plate_Location==7]
        Plate_Location <- 7
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==7]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==7]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==7]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==7]}
      }
      
      if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==8]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==8]
        Type<-METADat$Type[METADat$Plate_Location==8]
        Plate_Location <- 8
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==8]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==8]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==8]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==8]}
      }
      
      if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==9]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==9]
        Type<-METADat$Type[METADat$Plate_Location==9]
        Plate_Location <- 9
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==9]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==9]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==9]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==9]}
      }
      
      if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==10]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==10]
        Type<-METADat$Type[METADat$Plate_Location==10]
        Plate_Location <- 10
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==10]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==10]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==10]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==10]}
      }
      
      if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==11]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==11]
        Type<-METADat$Type[METADat$Plate_Location==11]
        Plate_Location <- 11
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==11]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==11]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==11]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==11]}
      }
      
      if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==12]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==12]
        Type<-METADat$Type[METADat$Plate_Location==12]
        Plate_Location <- 12
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==12]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==12]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==12]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==12]}
      }
      
      overall <- data.frame(filename = filename, Type= Type, well = well, row = row,
                            column = column, Dilution = Dilution, ID = ID, PA_Titer = PA_Titer,
                            Percent.SSC = Percent.SSC, Percent.FSC= Percent.FSC, Number_of_Cells=Number_of_Cells, Plate_Location = Plate_Location)
      return(overall)
    }
    
    z <- lapply(files,createoverall)
    O_Met_Dat <- rbind(z[[1]],z[[2]],z[[3]],z[[4]],z[[5]],z[[6]],z[[7]],z[[8]],z[[9]],z[[10]],z[[11]],z[[12]],z[[13]],z[[14]],z[[15]],z[[16]],z[[17]],z[[18]],z[[19]],z[[20]],z[[21]],z[[22]],z[[23]],z[[24]],z[[25]],z[[26]],z[[27]],z[[28]],z[[29]],z[[30]],z[[31]],z[[32]],z[[33]],z[[34]],z[[35]],z[[36]],z[[37]],z[[38]],z[[39]],z[[40]],z[[41]],z[[42]],z[[43]],z[[44]],z[[45]],z[[46]],z[[47]],z[[48]],z[[49]],z[[50]],z[[51]],z[[52]],z[[53]],z[[54]],z[[55]],z[[56]],z[[57]],z[[58]],z[[59]],z[[60]],z[[61]],z[[62]],z[[63]],z[[64]],z[[65]],z[[66]],z[[67]],z[[68]],z[[69]],z[[70]],z[[71]],z[[72]],z[[73]],z[[74]],z[[75]],z[[76]],z[[77]],z[[78]],z[[79]],z[[80]],z[[81]],z[[82]],z[[83]],z[[84]],z[[85]],z[[86]],z[[87]],z[[88]],z[[89]],z[[90]],z[[91]],z[[92]],z[[93]],z[[94]],z[[95]],z[[96]])
    
    O_Met_Dat$PA_Titer<-as.numeric(as.character(O_Met_Dat$PA_Titer))
    O_Met_Dat$Dilution<-as.numeric(as.character(O_Met_Dat$Dilution))
    O_Met_Dat$Act.Tit<-O_Met_Dat$PA_Titer/O_Met_Dat$Dilution
    O_Met_Dat$LN_Titer<-log(O_Met_Dat$Act.Tit)
    
    Negative<-mean(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    NegStDev<-sd(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    LOD<-3*NegStDev
    LOQ<-10*NegStDev
    O_Met_Dat$Adj.Per<-O_Met_Dat$Percent.SSC-Negative
    
    
    Calib<-lm(O_Met_Dat$LN_Titer[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ]~O_Met_Dat$Adj.Per[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ])    
    O_Met_Dat$FC_Titer[O_Met_Dat$Adj.Per>LOQ] <-exp(Calib$coefficients[1] + O_Met_Dat$Adj.Per[O_Met_Dat$Adj.Per>LOQ] * Calib$coefficients[2]) * O_Met_Dat$Dilution[O_Met_Dat$Adj.Per>LOQ]
    
                                                                                                   # plotting the calibration curve based on the "standard" blocks
    p<-ggplot(O_Met_Dat[O_Met_Dat$Type=="Standard",], aes(x=Percent.SSC, y=LN_Titer))
    p<-p + geom_point() + geom_smooth(method=lm, fullrange=TRUE)+xlim(0,100)+ylim(0,16)
    p <- print(p)
    curveInput <- p
  })
  
  persscdilInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    for (i in 1:length(Type)){
      if (Type[i] =="text/csv"){
        csvtext<- input$files$datapath[i]
      }
    }
    METADat <- read.csv(csvtext)
    O_Met_Dat <- 0
    
    createoverall <- function(filename) {
      flow <- read.FCS(file.path(filename))
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filename]

      
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      row<-gsub("[[:digit:]]", "", filename.split)
      column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))
      
      impDAT$FSC<-as.numeric(as.character(impDAT$FSC))
      
      modDAT<-impDAT[impDAT$FSC>=input$FSCgatelower & impDAT$FSC<=input$FSCgateupper,]
      Number_of_Cells<-length(impDAT$FSC)
      Percent.FSC<-length(modDAT$FSC)/length(impDAT$FSC)*100
      Percent.SSC<-length(which(as.numeric(as.character(modDAT$SSC))>input$ssclowlimit))/length(modDAT$FSC)*100

      #There might be a more succinct way of programming this section
      
      if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==1]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==1]
        Type<-METADat$Type[METADat$Plate_Location==1]
        Plate_Location <- 1
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==1]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==1]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==1]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==1]}
      }
      
      if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==2]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==2]
        Type<-METADat$Type[METADat$Plate_Location==2]
        Plate_Location <- 2
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==2]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==2]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==2]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==2]}
      }
      
      if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==3]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==3]
        Type<-METADat$Type[METADat$Plate_Location==3]
        Plate_Location <- 3
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==3]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==3]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==3]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==3]}
      }
      
      if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==4]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==4]
        Type<-METADat$Type[METADat$Plate_Location==4]
        Plate_Location <- 4
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==4]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==4]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==4]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==4]}
      }
      
      if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==5]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==5]
        Type<-METADat$Type[METADat$Plate_Location==5]
        Plate_Location <- 5
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==5]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==5]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==5]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==5]}
      }
      
      if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==6]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==6]
        Type<-METADat$Type[METADat$Plate_Location==6]
        Plate_Location <- 6
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==6]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==6]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==6]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==6]}
      }
      
      if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==7]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==7]
        Type<-METADat$Type[METADat$Plate_Location==7]
        Plate_Location <- 7
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==7]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==7]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==7]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==7]}
      }
      
      if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==8]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==8]
        Type<-METADat$Type[METADat$Plate_Location==8]
        Plate_Location <- 8
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==8]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==8]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==8]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==8]}
      }
      
      if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==9]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==9]
        Type<-METADat$Type[METADat$Plate_Location==9]
        Plate_Location <- 9
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==9]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==9]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==9]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==9]}
      }
      
      if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==10]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==10]
        Type<-METADat$Type[METADat$Plate_Location==10]
        Plate_Location <- 10
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==10]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==10]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==10]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==10]}
      }
      
      if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==11]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==11]
        Type<-METADat$Type[METADat$Plate_Location==11]
        Plate_Location <- 11
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==11]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==11]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==11]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==11]}
      }
      
      if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==12]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==12]
        Type<-METADat$Type[METADat$Plate_Location==12]
        Plate_Location <- 12
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==12]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==12]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==12]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==12]}
      }
      
      overall <- data.frame(filename = filename, Type= Type, well = well, row = row,
                            column = column, Dilution = Dilution, ID = ID, PA_Titer = PA_Titer,
                            Percent.SSC = Percent.SSC, Percent.FSC= Percent.FSC, Number_of_Cells=Number_of_Cells, Plate_Location = Plate_Location)
      return(overall)
    }
    
    z <- lapply(files,createoverall)
    O_Met_Dat <- rbind(z[[1]],z[[2]],z[[3]],z[[4]],z[[5]],z[[6]],z[[7]],z[[8]],z[[9]],z[[10]],z[[11]],z[[12]],z[[13]],z[[14]],z[[15]],z[[16]],z[[17]],z[[18]],z[[19]],z[[20]],z[[21]],z[[22]],z[[23]],z[[24]],z[[25]],z[[26]],z[[27]],z[[28]],z[[29]],z[[30]],z[[31]],z[[32]],z[[33]],z[[34]],z[[35]],z[[36]],z[[37]],z[[38]],z[[39]],z[[40]],z[[41]],z[[42]],z[[43]],z[[44]],z[[45]],z[[46]],z[[47]],z[[48]],z[[49]],z[[50]],z[[51]],z[[52]],z[[53]],z[[54]],z[[55]],z[[56]],z[[57]],z[[58]],z[[59]],z[[60]],z[[61]],z[[62]],z[[63]],z[[64]],z[[65]],z[[66]],z[[67]],z[[68]],z[[69]],z[[70]],z[[71]],z[[72]],z[[73]],z[[74]],z[[75]],z[[76]],z[[77]],z[[78]],z[[79]],z[[80]],z[[81]],z[[82]],z[[83]],z[[84]],z[[85]],z[[86]],z[[87]],z[[88]],z[[89]],z[[90]],z[[91]],z[[92]],z[[93]],z[[94]],z[[95]],z[[96]])
    
    
    O_Met_Dat$PA_Titer<-as.numeric(as.character(O_Met_Dat$PA_Titer))
    O_Met_Dat$Dilution<-as.numeric(as.character(O_Met_Dat$Dilution))
    O_Met_Dat$Act.Tit<-O_Met_Dat$PA_Titer/O_Met_Dat$Dilution
    O_Met_Dat$LN_Titer<-log(O_Met_Dat$Act.Tit)
    
    Negative<-mean(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    NegStDev<-sd(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    LOD<-3*NegStDev
    LOQ<-10*NegStDev
    O_Met_Dat$Adj.Per<-O_Met_Dat$Percent.SSC-Negative
    
    
    Calib<-lm(O_Met_Dat$LN_Titer[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ]~O_Met_Dat$Adj.Per[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ])    
    O_Met_Dat$FC_Titer[O_Met_Dat$Adj.Per>LOQ] <-exp(Calib$coefficients[1] + O_Met_Dat$Adj.Per[O_Met_Dat$Adj.Per>LOQ] * Calib$coefficients[2]) * O_Met_Dat$Dilution[O_Met_Dat$Adj.Per>LOQ]
    
                                                                                                   #plotting percent ssc as a function of dilution for each of the blocks marked sample
    p<-ggplot(O_Met_Dat[O_Met_Dat$Type=="Sample",], aes(x=Dilution, y=Percent.SSC))
    p<-p + geom_point() + geom_smooth(method=lm)
    p<-p + facet_wrap(~ID)
    print(p)
    
    persscdilInput <- p
  })
  
  
  fscsscInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
                                                                                                   #functionally similar to O_Met_Dat in other portions of the code, completefscssc holds information of all ssc and fsc VALUES for all 96 wells
    completefscssc <- 0
    
    for (filepath in files) {                                                                      
      flow <- read.FCS(file.path(filepath))                                                       
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filepath]                              
      
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      total <- cbind(rep(well,length(impDAT[[1]])),impDAT)
      if (completefscssc==0) {completefscssc <- total}                                             
      else (completefscssc <- rbind(completefscssc,total))
    }
    colnames(completefscssc) <- c('Well','SSC','FSC')
 
                                                                                                   #plots ssc values compared to fsc values for cells in all 96 wells. takes a while to load
    p <- ggplot(completefscssc, aes(x=FSC, y=SSC))+
    geom_point(alpha=0.05)+ 
    facet_wrap(~Well, nrow=8)+
    scale_y_log10()+ 
    scale_x_log10()
    print(p)
    
    fscsscInput <- p
  })
  
  negsdlodloqInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    for (i in 1:length(Type)){
      if (Type[i] =="text/csv"){
        csvtext<- input$files$datapath[i]
      }
    }
    METADat <- read.csv(csvtext)
    O_Met_Dat <- 0
    
    createoverall <- function(filepath) {
      flow <- read.FCS(file.path(filepath))
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filepath]
      
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      row<-gsub("[[:digit:]]", "", filename.split)
      column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))

      impDAT$FSC<-as.numeric(as.character(impDAT$FSC))
      
      modDAT<-impDAT[impDAT$FSC>=input$FSCgatelower & impDAT$FSC<=input$FSCgateupper,]
      Number_of_Cells<-length(impDAT$FSC)
      Percent.FSC<-length(modDAT$FSC)/length(impDAT$FSC)*100
      Percent.SSC<-length(which(as.numeric(as.character(modDAT$SSC))>input$ssclowlimit))/length(modDAT$FSC)*100
      
      if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==1]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==1]
        Type<-METADat$Type[METADat$Plate_Location==1]
        Plate_Location <- 1
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==1]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==1]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==1]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==1]}
      }
      
      if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==2]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==2]
        Type<-METADat$Type[METADat$Plate_Location==2]
        Plate_Location <- 2
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==2]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==2]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==2]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==2]}
      }
      
      if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==3]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==3]
        Type<-METADat$Type[METADat$Plate_Location==3]
        Plate_Location <- 3
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==3]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==3]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==3]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==3]}
      }
      
      if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==4]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==4]
        Type<-METADat$Type[METADat$Plate_Location==4]
        Plate_Location <- 4
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==4]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==4]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==4]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==4]}
      }
      
      if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==5]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==5]
        Type<-METADat$Type[METADat$Plate_Location==5]
        Plate_Location <- 5
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==5]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==5]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==5]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==5]}
      }
      
      if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
        ID<- METADat$ID[METADat$Plate_Location==6]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==6]
        Type<-METADat$Type[METADat$Plate_Location==6]
        Plate_Location <- 6
        if (row=="A"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==6]}
        if (row=="B"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==6]}
        if (row=="C"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==6]}
        if (row=="D"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==6]}
      }
      
      if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==7]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==7]
        Type<-METADat$Type[METADat$Plate_Location==7]
        Plate_Location <- 7
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==7]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==7]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==7]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==7]}
      }
      
      if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==8]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==8]
        Type<-METADat$Type[METADat$Plate_Location==8]
        Plate_Location <- 8
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==8]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==8]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==8]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==8]}
      }
      
      if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==9]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==9]
        Type<-METADat$Type[METADat$Plate_Location==9]
        Plate_Location <- 9
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==9]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==9]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==9]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==9]}
      }
      
      if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==10]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==10]
        Type<-METADat$Type[METADat$Plate_Location==10]
        Plate_Location <- 10
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==10]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==10]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==10]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==10]}
      }
      
      if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==11]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==11]
        Type<-METADat$Type[METADat$Plate_Location==11]
        Plate_Location <- 11
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==11]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==11]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==11]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==11]}
      }
      
      if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
        ID<- METADat$ID[METADat$Plate_Location==12]
        PA_Titer<-METADat$PA_Titer[METADat$Plate_Location==12]
        Type<-METADat$Type[METADat$Plate_Location==12]
        Plate_Location <- 12
        if (row=="E"){Dilution<-METADat$Dilution_1[METADat$Plate_Location==12]}
        if (row=="F"){Dilution<-METADat$Dilution_2[METADat$Plate_Location==12]}
        if (row=="G"){Dilution<-METADat$Dilution_3[METADat$Plate_Location==12]}
        if (row=="H"){Dilution<-METADat$Dilution_4[METADat$Plate_Location==12]}
      }
      
      overall <- data.frame(filename = filename, Type= Type, well = well, row = row,
                            column = column, Dilution = Dilution, ID = ID, PA_Titer = PA_Titer,
                            Percent.SSC = Percent.SSC, Percent.FSC= Percent.FSC, Number_of_Cells=Number_of_Cells, Plate_Location = Plate_Location)
      return(overall)
    }
    
    z <- lapply(files,createoverall)
    O_Met_Dat <- rbind(z[[1]],z[[2]],z[[3]],z[[4]],z[[5]],z[[6]],z[[7]],z[[8]],z[[9]],z[[10]],z[[11]],z[[12]],z[[13]],z[[14]],z[[15]],z[[16]],z[[17]],z[[18]],z[[19]],z[[20]],z[[21]],z[[22]],z[[23]],z[[24]],z[[25]],z[[26]],z[[27]],z[[28]],z[[29]],z[[30]],z[[31]],z[[32]],z[[33]],z[[34]],z[[35]],z[[36]],z[[37]],z[[38]],z[[39]],z[[40]],z[[41]],z[[42]],z[[43]],z[[44]],z[[45]],z[[46]],z[[47]],z[[48]],z[[49]],z[[50]],z[[51]],z[[52]],z[[53]],z[[54]],z[[55]],z[[56]],z[[57]],z[[58]],z[[59]],z[[60]],z[[61]],z[[62]],z[[63]],z[[64]],z[[65]],z[[66]],z[[67]],z[[68]],z[[69]],z[[70]],z[[71]],z[[72]],z[[73]],z[[74]],z[[75]],z[[76]],z[[77]],z[[78]],z[[79]],z[[80]],z[[81]],z[[82]],z[[83]],z[[84]],z[[85]],z[[86]],z[[87]],z[[88]],z[[89]],z[[90]],z[[91]],z[[92]],z[[93]],z[[94]],z[[95]],z[[96]])
    
    O_Met_Dat$PA_Titer<-as.numeric(as.character(O_Met_Dat$PA_Titer))
    O_Met_Dat$Dilution<-as.numeric(as.character(O_Met_Dat$Dilution))
    O_Met_Dat$Act.Tit<-O_Met_Dat$PA_Titer/O_Met_Dat$Dilution
    O_Met_Dat$LN_Titer<-log(O_Met_Dat$Act.Tit)
    
    Negative<-mean(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    NegStDev<-sd(O_Met_Dat$Percent.SSC[O_Met_Dat$Type=="Blank"])
    LOD<-3*NegStDev
    LOQ<-10*NegStDev
    O_Met_Dat$Adj.Per<-O_Met_Dat$Percent.SSC-Negative
    
    
    Calib<-lm(O_Met_Dat$LN_Titer[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ]~O_Met_Dat$Adj.Per[O_Met_Dat$Type=="Standard" & O_Met_Dat$Adj.Per>LOQ])    
    O_Met_Dat$FC_Titer[O_Met_Dat$Adj.Per>LOQ] <-exp(Calib$coefficients[1] + O_Met_Dat$Adj.Per[O_Met_Dat$Adj.Per>LOQ] * Calib$coefficients[2]) * O_Met_Dat$Dilution[O_Met_Dat$Adj.Per>LOQ]
    
                                                                                                   # displays the negative standard deviation, limit of detection, and limit of quantification
    negsdlodloq <- cbind(as.data.frame(c("Negative Standard Deviation","Limit of Detection","Limit of Quantification")), as.data.frame(c(NegStDev,LOD,LOQ)))
    colnames(negsdlodloq) <- c("Data","Values")
    
    negsdlodloqlInput <- negsdlodloq
  })
  
  fscsschistInput <- reactive ({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    fscsscvalues <- 0                                             
    
    for (filepath in files) {                                      
      flow <- read.FCS(file.path(filepath))                         
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filepath] 
      
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      total <- cbind(rep(well,length(impDAT[[1]])),impDAT)
      if (fscsscvalues==0) {fscsscvalues <- total}            
      else (fscsscvalues <- rbind(fscsscvalues,total))
    }
    
    colnames(fscsscvalues) <- c('Well','SSC','FSC')
    types1 <- rep('SSC',length(fscsscvalues$SSC))
    types2 <- rep('FSC',length(fscsscvalues$FSC))
    typestotal <- c(types1,types2)
    totaldata <- cbind(as.data.frame(typestotal),as.data.frame(c(fscsscvalues$SSC,fscsscvalues$FSC)))
    colnames(totaldata) <- c('Type','Value')
  
                                                                                                   #plots 2 histograms showing the total fsc and ssc values
    g <- ggplot(data=totaldata,aes(x=Value))+geom_histogram(bins=input$histbins)+facet_grid(Type~.)
    print(g)
    fscsschistInput <- g
  })
  
  sscblankplotInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    for (i in 1:length(Type)){
      if (Type[i] =="text/csv"){
        csvtext<- input$files$datapath[i]
      }
    }
    METADat <- read.csv(csvtext)
    locssc <- 0
    
    for (filepath in files) {
      flow <- read.FCS(file.path(filepath))
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
      
      filename <- name.and.files$name[name.and.files$files==filepath]
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      row<-gsub("[[:digit:]]", "", filename.split)
      column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))

      if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 1
        Type<-METADat$Type[METADat$Plate_Location==1]
      }
      if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 2
        Type<-METADat$Type[METADat$Plate_Location==2]
      }
      if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 3
        Type<-METADat$Type[METADat$Plate_Location==3]
      }
      if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 4
        Type<-METADat$Type[METADat$Plate_Location==4]
      }
      if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 5
        Type<-METADat$Type[METADat$Plate_Location==5]
      }
      if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 6
        Type<-METADat$Type[METADat$Plate_Location==6]
      }
      if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 7
        Type<-METADat$Type[METADat$Plate_Location==7]
      }
      if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 8
        Type<-METADat$Type[METADat$Plate_Location==8]
      }
      if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 9
        Type<-METADat$Type[METADat$Plate_Location==9]
      }
      if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 10
        Type<-METADat$Type[METADat$Plate_Location==10]
      }
      if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 11
        Type<-METADat$Type[METADat$Plate_Location==11]
      }
      if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 12
        Type<-METADat$Type[METADat$Plate_Location==12]
      }
      
      location.value <- cbind(rep(Plate_Location,length(impDAT$SSC)),impDAT$SSC)
      if (Type=="Blank") {
        if (locssc==0){locssc<-location.value}
        else {locssc<-rbind(locssc,location.value)}
      }

    }
    locssc <- as.data.frame(locssc)
    colnames(locssc)<- c("Location","Value")

    graphs <- unique(locssc$Location)
    
                                                                                                   #plots multiple distribution curves based on the number of blank blocks there are, with a vertical line indicating where the lower ssc limit is set to
    p <- ggplot(locssc,aes(x=Value,color=as.factor(Location))) +
      geom_density(aes(color=as.factor(Location)))+
      geom_vline(xintercept=input$ssclowlimit,color="red",size = 1.3, alpha = 0.5)

    sscblankplotInput <- print(p)
  })
  
  sscblanktableInput <- reactive({
    Type <- input$files$type
    
    files <- vector()
    name <- vector()
    
    for (i in 1:length(Type)) {
      if (Type[i] !="text/csv"){
        files <- c(files,input$files$datapath[i])
        name <- c(name,input$files$name[i])
      }
    }
    name.and.files <- cbind(as.data.frame(name),as.data.frame(files))
    
    for (i in 1:length(Type)){
      if (Type[i] =="text/csv"){
        csvtext<- input$files$datapath[i]
      }
    }
    METADat <- read.csv(csvtext)
    locssc <- 0
    
    for (filepath in files) {
      flow <- read.FCS(file.path(filepath))
      impDAT <- data.frame(flow@exprs)[ , c('SSC.A', 'FSC.A')]
      colnames(impDAT) <- c('SSC','FSC')
  
      filename <- name.and.files$name[name.and.files$files==filepath]
      filename.split<- gsub('.fcs','',filename)
      well<-filename.split
      row<-gsub("[[:digit:]]", "", filename.split)
      column<-as.numeric(gsub("[^[:digit:]]", "", filename.split))
      
      if ((column == 1 | column == 2) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 1
        Type<-METADat$Type[METADat$Plate_Location==1]
      }
      if ((column == 3 | column == 4) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 2
        Type<-METADat$Type[METADat$Plate_Location==2]
      }
      if ((column == 5 | column == 6) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 3
        Type<-METADat$Type[METADat$Plate_Location==3]
      }
      if ((column == 7 | column == 8) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 4
        Type<-METADat$Type[METADat$Plate_Location==4]
      }
      if ((column == 9 | column == 10) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 5
        Type<-METADat$Type[METADat$Plate_Location==5]
      }
      if ((column == 11 | column == 12) & (row=="A" | row=="B" | row=="C" | row=="D")){
        Plate_Location <- 6
        Type<-METADat$Type[METADat$Plate_Location==6]
      }
      if ((column == 1 | column == 2) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 7
        Type<-METADat$Type[METADat$Plate_Location==7]
      }
      if ((column == 3 | column == 4) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 8
        Type<-METADat$Type[METADat$Plate_Location==8]
      }
      if ((column == 5 | column == 6) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 9
        Type<-METADat$Type[METADat$Plate_Location==9]
      }
      if ((column == 7 | column == 8) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 10
        Type<-METADat$Type[METADat$Plate_Location==10]
      }
      if ((column == 9 | column == 10) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 11
        Type<-METADat$Type[METADat$Plate_Location==11]
      }
      if ((column == 11 | column == 12) & (row=="E" | row=="F" | row=="G" | row=="H")){
        Plate_Location <- 12
        Type<-METADat$Type[METADat$Plate_Location==12]
      }
      
                                                                                                   #for the blank blocks, makes a dataframe with the block number next to each ssc value
      location.value <- cbind(rep(Plate_Location,length(impDAT$SSC)),impDAT$SSC)
      if (Type=="Blank") {
        if (locssc==0){locssc<-location.value}
        else {locssc<-rbind(locssc,location.value)}
      }
    }
    
    locssc <- as.data.frame(locssc)
    colnames(locssc)<- c("Location","Value")
    blocks <- unique(locssc$Location)
    
                                                                                             
    table <- matrix()
    
                                                                                                   #for 1 to the number of blank blocks,
    for (i in 1:length(blocks)) {
      
                                                                                                   #calculates the total number of ssc values for the current block
      blocklength <- length(locssc$Value[locssc$Location==blocks[i]])
      
                                                                                                   #calculates the percent of ssc values in the block that are above the ssc lower limit
      percent <- (length(which(locssc$Value[locssc$Location==blocks[i]] > input$ssclowlimit))/length(locssc$Value[locssc$Location==blocks[i]]))*100
      table[i] <- percent
      
    }
    
                                                                                                   #makes a flat table with the blocks that are blank over their corresponding percent ssc value
    table <- as.data.frame(t(table))
    row.names(table) <- "Percent"
    colnames(table) <- as.character(blocks)
    
    sscblanktableInput <- table
  })
############################################################################### Outputs ###########################################################################################                         
   
  # linking the dataInput portion to the ID "dataOutput" by calling it in the dataOutput section of the output section of the server so it shows up when the part on the main frame displaying the dataOutput table 
  output$dataOutput <- renderTable({                      
     if (is.null(input$files))
       return(NULL)
     dataInput()
   })
   
   output$curveOutput <- renderPlot({
     if (is.null(input$files))
       return(NULL)
     curveInput()
   })
   
   output$persscdilOutput <- renderPlot({
     if (is.null(input$files))
       return(NULL)
     persscdilInput()
   })
   
   output$fscsscOutput <- renderPlot({
     if (is.null(input$files))
       return(NULL)
     fscsscInput()
   },
   width = 1400,
   height=1000)
   
   output$negsdlodloqOutput <- renderTable({
     if (is.null(input$files))
       return(NULL)
     negsdlodloqInput()
   })
   
   output$fscsschistOutput <- renderPlot({
     if (is.null(input$files))
       return(NULL)
     fscsschistInput()
   })
   
   output$sscblankplotOutput <- renderPlot({
     if (is.null(input$files))
       return(NULL)
     sscblankplotInput()
   })
   output$sscblanktableOutput <- renderTable({
     if(is.null(input$files))
       return(NULL)
     sscblanktableInput()
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

