#' @title The smallest number.
#' @name minR
#' @description Find the smallest number in a vector
#' @useDynLib StatComp20041
#' @examples
#' \dontrun{
#' x=c(rep(1,100),rep(2,100),rep(0,100))
#' minR(x)
#' }
#' @export
minR<-function(x){
  num<-length(x)
  if(num==1){return(x)
  }else{
    x_min<-x[1]
    for(i in 2:length(x)){
      if(x_min>x[i]){x_min<-x[i]}
    }
    return(x_min)
  }
}
#' @title Choose sort method
#' @name choose_sort
#' @description Use choose sort method to rank a vector
#' @useDynLib StatComp20041
#' @examples
#' \dontrun{
#' x=c(rep(1,100),rep(2,100),rep(0,100))
#' choose_sort(x)
#' }
#' @export
choose_sort<-function(x){
  u<-x
  s<-c()
  while(length(u)>0){
    u_min<-minR(u)
    s<-c(s,u_min)
    u<-u[-which(u==u_min)[1]]
  }
  return(s)
}
#' @title Insert sort method
#' @name insert_sort
#' @description Use insert sort method to rank a vector
#' @useDynLib StatComp20041
#' @examples
#' \dontrun{
#' x=c(rep(1,100),rep(2,100),rep(0,100))
#' insert_sort method(x)
#' }
#' @export
insert_sort<-function(x){
  if(length(x)==1){return(x)}
  u<-x[-length(x)]
  s<-x[length(x)]
  while(length(u)>0){
    n<-length(u)
    s1<-s[s-u[n]<0]
    s2<-s[s-u[n]>=0]
    s<-c(s1,u[n],s2)#保持s的顺序
    u<-u[-n]
  }
  return(s)
}
#' @title Bubble sort method
#' @name bubble_sort
#' @description Use bubble sort method to rank a vector
#' @useDynLib StatComp20041
#' @examples
#' \dontrun{
#' x=c(rep(1,100),rep(2,100),rep(0,100))
#' bubble_sort(x)
#' }
#' @export
bubble_sort<-function(x){
  num<-length(x)
  if(num==1){return(x)
  }else{
    finished<-FALSE
    while(!finished){
      m<-x
      for(i in 1:(num-1)){
        s<-x
        if(x[i]>x[i+1]){
          s[i]<-x[i+1]
          s[i+1]<-x[i]}
        x<-s
      }
      finished<-sum(s==m)==num#判断是否不存在需要交换的元素
    }
    return(x)
  } 
}
#' @title Quick sort method
#' @name quick_sort
#' @description Use quick sort method to rank a vector
#' @useDynLib StatComp20041
#' @examples
#' \dontrun{
#' x=c(rep(1,100),rep(2,100),rep(0,100))
#' quick_sort(x)
#' }
#' @export
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#递归
}
