get_last <- function(data){
myid.uni <- unique(data$user_id)
a<-length(myid.uni)

last <- c()

for (i in 1:a) {
  temp<-subset(data, user_id==myid.uni[i])
  if (dim(temp)[1] > 1) {
    last.temp<-temp[dim(temp)[1],]
  }
  else {
    last.temp<-temp
  }
  last<-rbind(last, last.temp)
}
last
}