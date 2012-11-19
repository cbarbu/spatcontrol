library(testthat)
maps<-read.csv("maps_hunter.csv")
blockMap<-read.csv("Hunter_Points_Casas-blocks.csv.R")

blocksLinesInBlockMap<-match(maps$unicode_gps,blockmap$unicode)
maps$block_num<-blockmap$polygon[blocksLinesInBlockMap]

# following should be a map of colored city-blocks
with(maps,plot(X,Y,col=block_num,asp=1))

# following should not return an error message 
expect_equal(length(which(maps$block_num==-1)),0)

write.csv(maps,"maps_hunter_blocks.csv",row.names=FALSE)
