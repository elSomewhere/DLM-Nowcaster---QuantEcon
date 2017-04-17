dbDisconnect(db)
print("disconnected")


rm(list=ls())
setwd("F:/FI-FX/RESEARCH_INTERN/Esteban/nowcasting/Nowcast_November")

library(openxlsx)
library(Rblpapi)
library(zoo)
library(RSQLite)
blpConnect()

source("NC_FrameworkUtils.R")
source("DLMMaster.R")
source("NEWS_DLMKalmanFilter_bySQL.R")
source("NEWS_DLMKalmanFilter_byList.R")
source("NEWS_computeIterativeNews_7_test.R")

#####################
#===================#
#####################
rebuildDataBase <- TRUE
recycleParameters <- TRUE
computeNews <- TRUE
rebuildAllNews <- FALSE
country <- "US"
dataSource <- "BBG"
maxForecast <- 3
maxNewInfoHorizon <- 3


############################
#prepare database of models#
############################
dbName <- paste0(country,"_",dataSource,"_Output.sqlite")
if(rebuildDataBase==TRUE){
  print("re-initializing database - notice that all old model records get deleted")
  file.remove(dbName)
  drive <- SQLite()
  db <- dbConnect(drive, dbname=dbName)
  dbGetQuery(conn = db, 
             "CREATE TABLE IF NOT EXISTS models 
             (_id INTEGER PRIMARY KEY AUTOINCREMENT, 
             date TEXT UNIQUE,
             model BLOB)")
  dbGetQuery(conn = db, 
             "CREATE TABLE IF NOT EXISTS news 
             (_id INTEGER PRIMARY KEY AUTOINCREMENT, 
             date TEXT UNIQUE,
             newStuff BLOB,
             FOREIGN KEY(date) REFERENCES models(date)
             )")
}else{
  print("connecting to existing database and appending output")
  drive <- SQLite()
  db <- dbConnect(drive, dbname=dbName)
}

queryNumRows <- "SELECT _id FROM models"
queryNumRows <- dbGetQuery(db,queryNumRows)
queryNumRows <- length(queryNumRows$`_id`)
print(paste0("there are currently ", queryNumRows, " rowentries in the data base"))



##############
#collect data#
##############
specs_raw <- read.xlsx("SPECS_BBG.XLSX", sheet = country, startRow = 1, colNames = TRUE,rowNames = FALSE)
specs <- lapply(seq(dim(specs_raw)[1]),function(x){
  res <- list(specs_raw[x,2],specs_raw[x,3],specs_raw[x,4],specs_raw[x,5],specs_raw[x,6],specs_raw[x,7],specs_raw[x,8],specs_raw[x,9],specs_raw[x,10])
  names(res) <- colnames(specs_raw)[-1]
  return(res)
})
names(specs) <- specs_raw[,1]

start <- as.Date("1985-01-31")
end <- Sys.Date()+45
today <- Sys.Date()

retrieval <- lapply(specs,function(x){
  print(x["BBG"])
  naked <- getBBG(x["BBG"][[1]],start,end,x["Freq"][[1]])[,-2]
  eval <- x["Transform"]
  if(eval == 1){
    transformed <- zoo2Function(zoo2Function((naked),takeLog),takeDiff,1)
  }else if(eval==2){
    transformed <- zoo2Function((naked),takeDiff,1)
  }else if(eval==3){
    transformed <- zoo2Function(zoo2Function(zoo2Function(zoo2Function(naked,cumSeries),raiseToPos),takeLog),takeDiff,1)
  }else if(eval==0){
    transformed <- naked
  }
  naked <- zoo(naked,order.by=as.yearmon(as.Date(index(naked),frac=1)))
  transformed <- zoo(transformed,order.by=as.yearmon(as.Date(index(transformed),frac=1)))
  res <- list(naked,transformed)
  names(res) <- c("naked","transformed")
  return(res)
})
naked <- lapply(retrieval,function(x){
  res <- x[[1]]
  return(res)
})
transformed <- lapply(retrieval,function(x){
  res <- x[[2]]
  return(res)
})
data_transformed <- do.call(cbind,transformed)
data_naked <- do.call(cbind,naked)



############
#clean data#
############
data <- easyTrimmer_soft(data_transformed,1)
dates <- index(data)
names <- colnames(data)
data_raw <- as.matrix(data)
colnames(data_raw) <- NULL
rownames(data_raw) <- NULL


############
#write data#
############
write.csv(data_raw,"USDATA.csv")
wb <- createWorkbook()
addWorksheet(wb, sheetName = "USDATA", gridLines = FALSE)
writeData(wb, sheet = "USDATA", x = cbind(index(data),data),colNames=TRUE,rowNames=TRUE)
saveWorkbook(wb, "USDATA.xlsx", overwrite = TRUE)

#####################
#===================#
#####################
index_idio_m <- t(matrix(which(specs_raw$Freq=="M")))
index_idio_q <- t(matrix(which(specs_raw$Freq=="Q")))

#442 rolling
#423 recursive

r <- 4
p <- 2
g <- 3
fch <- 12
bch <- 12


###############
#Compute Model#
###############
if(queryNumRows==0||recycleParameters==FALSE){
  print("re-estimating all the parameters")
  testDLFM <- DLMCockpit(data_raw,r,p,g,index_idio_m,index_idio_q,preInitSteps=10,fch,bch,max_iter=25,MLE=FALSE)
}else{
  print("using parameters of last model solved for intialization")
  query <- "SELECT * FROM models WHERE _id=(SELECT MAX (_id) FROM models)"
  lastModel <- dbGetQuery(db,query)
  lastModel <- unserialize(lastModel$model[[1]])
  lastModel <- lastModel$model$model_opt
  testDLFM <- DLMShortcut(data_raw,r,p,g,index_idio_m,index_idio_q,lastModel,max_iter=25,bch=bch,fch=fch,MLE=FALSE,reOptimize=FALSE)
}




##################
#transform output#
##################

dates_fc <- c(dates[1]-rev((seq(bch)/12)),dates,dates[length(dates)]+((seq(fch)/12)))
output_raw <- testDLFM$outPutSeries
output_raw_zoo <- zoo(testDLFM$outPutSeries,order.by=dates_fc)
colnames(output_raw_zoo) <- specs_raw[,1]    


print(data[as.yearmon(Sys.Date()),])
print(output_raw_zoo[as.yearmon(Sys.Date()),])

model <- list(testDLFM,names,data,output_raw_zoo)
names(model) <- c("model","names","rawInput_zoo","rawOutput_zoo")


#####################
#Write output charts#
#####################

writeOutput(testDLFM,names,dates_fc)

model <- list(testDLFM,names,dates_fc)
names(model) <- c("model","names","dates")
print("finished writing output")













#######################
#Save data to database#
#######################
queryDates <- "SELECT date FROM models"
queryDates <- dbGetQuery(db,queryDates)
queryDates <- as.matrix(queryDates)

if(length(queryDates)==0){
  print("created the first entry into the db")
}else{
  if(queryDates[dim(queryDates)[1],]==today){
    print("entry for today already exists, replacing that entry with this one")
  }else{
    print("appended a new day into the database")
  }
}

forOutput <- serialize(model,NULL)
toSQL <- data.frame(d=as.character(today),g=I(list(forOutput)))
dbGetPreparedQuery(db, "INSERT OR REPLACE INTO models (date,model) VALUES (:d,:g)", bind.data=toSQL)

queryNumRowsNew <- "SELECT _id FROM models"
queryNumRowsNew <- dbGetQuery(db,queryNumRowsNew)
queryNumRowsNew <- length(queryNumRowsNew$`_id`)
print(paste0("there are now currently ", queryNumRowsNew, " rowentries in the data base"))


##############
#compute News#
##############
numModelentries <- dim(dbGetQuery(conn = db,"SELECT * FROM models"))[1]
numNewsEntries <- dim(dbGetQuery(conn = db,"SELECT * FROM news"))[1]
modelIDs <- dbGetQuery(db,"SELECT _id FROM models")
modelIDs <- modelIDs$`_id`
modelDates <- dbGetQuery(conn = db,"SELECT date FROM models")
modelDates <- modelDates$date

if(computeNews==TRUE&&numModelentries>1){
  if(rebuildAllNews==TRUE){
    dbGetQuery(conn = db,"DELETE FROM news")
    dbGetQuery(conn = db,"DELETE FROM SQLITE_SEQUENCE where name='news'")
    
    for(i in 2:numModelentries){
      print(i)
      id_post <- modelIDs[i]
      id_pre <- modelIDs[i-1]
      pre <- paste0("SELECT * FROM models WHERE _id = ",id_pre)
      pre <- dbGetQuery(db,pre)
      pre <- unserialize(pre$model[[1]])
      post <- paste0("SELECT * FROM models WHERE _id = ",id_post)
      post <- dbGetQuery(db,post)
      postDate <- post$date
      post <- unserialize(post$model[[1]])
      
      ##############
      #compute news#
      ##############
      news <- computeNews_byMem(pre,post,maxForecast=maxForecast,maxNewInfoHorizon=maxNewInfoHorizon)
      #news <- computeNews_bySQL(pre,post,maxForecast=3,maxNewInfoHorizon=3)
      
      #######################
      #Save data to database#
      #######################
      forOutput <- serialize(news,NULL)
      toSQL <- data.frame(d=as.character(postDate),g=I(list(forOutput)))
      dbGetPreparedQuery(db, "INSERT OR REPLACE INTO news (date,newStuff) VALUES (:d,:g)", bind.data=toSQL)
    }
  }else{
    #just append the newest
    i <- numModelentries
    print(i)
    id_post <- modelIDs[i]
    id_pre <- modelIDs[i-1]
    pre <- paste0("SELECT * FROM models WHERE _id = ",id_pre)
    pre <- dbGetQuery(db,pre)
    pre <- unserialize(pre$model[[1]])
    post <- paste0("SELECT * FROM models WHERE _id = ",id_post)
    post <- dbGetQuery(db,post)
    postDate <- post$date
    post <- unserialize(post$model[[1]])
    
    ##############
    #compute news#
    ##############
    news <- computeNews_byMem(pre,post,maxForecast=maxForecast,maxNewInfoHorizon=maxNewInfoHorizon)
    #news <- computeNews_bySQL(pre,post,maxForecast=3,maxNewInfoHorizon=3)
    
    #######################
    #Save data to database#
    #######################
    print("appending news to database")
    forOutput <- serialize(news,NULL)
    toSQL <- data.frame(d=as.character(postDate),g=I(list(forOutput)))
    dbGetPreparedQuery(db, "INSERT OR REPLACE INTO news (date,newStuff) VALUES (:d,:g)", bind.data=toSQL)
  }
}else{
  print("cannot comupte news yet, need more models in database (at least 2...")
}
}

############
#Disconnect#
############
dbDisconnect(db)
print("disconnected")



