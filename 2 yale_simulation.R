noi <- j <- 0.3 #P=0.6
n1<-n2<-100;n3<-3
train <- train.noi <- array(0,c(n1,n2,n3))
for(i in 1:n3){
  train[,,i] <- read.bmp(sprintf("Yaledata/0%d.bmp", i))
  train[,,i] <- (train[,,i]-min(train[,,i]))/(max(train[,,i])-min(train[,,i]))
  train.noi[,,i] <-  train[,,i]
  train.noi[,,i][sample(which(train.noi[,,i]>0),100*100*j)]<-0
  train.noi[,,i][sample(which(train.noi[,,i]>0),100*100*j)]<-1
}
#导出噪声图片
for(i in 1:n3){
  writePNG(train.noi[,,i],sprintf("noise/0%d.png", i))
}
#mine
pro <- array(0,c(n1,n2,n3))
pro <- Nee(train.noi,10/sqrt(100),18/sqrt(100),2,50,50,0.5)
for(i in 1:n3){
  writePNG(pro[,,i],sprintf("result/0%d.png", i))
}
