setwd("C:/Users/Esteban/Dropbox/barcelona/nowcaster")
library(RSQLite)
model1 <- readRDS("model1.rds")
model2 <- readRDS("model2.rds")
model3 <- readRDS("model3.rds")
out1 <- serialize(model1,NULL)
out2 <- serialize(model2,NULL)
out3 <- serialize(model3,NULL)

db <- dbConnect(SQLite(), dbname="Test2.sqlite")
dbGetQuery(conn = db, 
           "CREATE TABLE IF NOT EXISTS models 
            (_id INTEGER PRIMARY KEY AUTOINCREMENT, 
            date TEXT,
            model BLOB)")

test1 <- data.frame(d="1980-01-01",g=I(list(out1)))
test2 <- data.frame(d="1982-01-01",g=I(list(out2)))
test3 <- data.frame(d="1983-01-01",g=I(list(out3)))
dbGetPreparedQuery(db, "INSERT INTO models (date,model) VALUES (:d,:g)", bind.data=test1)
dbGetPreparedQuery(db, "INSERT INTO models (date,model) VALUES (:d,:g)", bind.data=test2)
dbGetPreparedQuery(db, "INSERT INTO models (date,model) VALUES (:d,:g)", bind.data=test3)

dbListTables(db)

p1 <- dbGetQuery( db,'select * from models' )
p2 <- dbGetQuery(db, "SELECT * FROM models WHERE Date(date) >= Date('1982-01-01')")
testin1 <- unserialize(p1$model[[3]])
testin2 <- unserialize(p1$model[[1]])