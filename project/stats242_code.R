library(mhsmm)

get.stock.data.weekly = function(stock.name, data) {
  data.list = list();
  for (s in stock.name) {
    t = data[data$TICKER==s, ];
    t$IND = (1:length(t$DATE))
    t = t[t$IND %% 5 == 1, ];
    t$RET = c(0, diff(log(t$PRC)));
    t$IND = (1:length(t$IND));
    data.list[[length(data.list) + 1]] = t;
  }  
  return(data.list);
}


get.stock.data = function(stock.name, data) {
  data.list = list();
  for (s in stock.name) {
    t = data[data$TICKER==s, ];
    t$RET = c(0, diff(log(t$PRC)));
    data.list[[length(data.list) + 1]] = t;
  }  
  return(data.list);
}

pick.num.clusters = function(stock) {
  ret = stock$RET;
  ret = ret[2:length(ret)];
  name = stock$TICKER[1];
  x = rep(0,7)
  for (k in 1:7){
    km = kmeans(ret, centers = k);
    x[k] = km$betweenss / km$totss;	
  }
  plot(x,type='o',main = name, xlab="Number of cluster",
       ylab="Percentage of Variance Explained");
}

get.initial.gaussian = function(ret, num.clusters) {
  km = kmeans(ret, centers = num.clusters);
  variance = c();
  for (i in 1:num.clusters) {
    t = ret[km$cluster == i]
    variance[i] = var(t);
  }
  return(list(mu=c(km$centers), sigma=variance));
}

is.good.state = function(stock, num.clusters, cur.ind, num.look.back) {
  stock = stock[(cur.ind-num.look.back+1):cur.ind, ];
  ret = stock$RET;
  ret = ret[2:length(ret)];
  init0 = rep(1/num.clusters, num.clusters);
  P0 = matrix(1/num.clusters, nrow=num.clusters, ncol=num.clusters);
  b0 = get.initial.gaussian(ret, num.clusters);
  startval = hmmspec(init=init0, trans=P0, parms.emission=b0, dens.emission=dnorm.hsmm);
  h = hmmfit(list(x=ret, N=length(ret)), startval, mstep=mstep.norm, tol=1e-8, maxit=10000);
  cur.state = h$yhat[length(h$yhat)];
  mu = h$model$parms.emission$mu;
  sigma = h$model$parms.emission$sigma;
  if (mu[cur.state] - sigma[cur.state] > 0) {
    return (TRUE);
  } else {
    return (FALSE);
  }
}

trade.buy.and.hold = function(stock, start.ind) {
  end.ind = length(stock$RET)
  ret.t = c();
  hold = FALSE;
  for (i in start.ind:end.ind) {
    if (hold == TRUE) {
      ret.t = c(ret.t, 1 + stock$RET[i])
    } else {
      ret.t = c(ret.t, 1);
    }
    hold = TRUE;
  }
  return(ret.t)
}

trade.hmm = function(stock, look.back, num.clusters) {
  start.ind = look.back;
  end.ind = length(stock$RET);
  ret.t = c();
  prev.state = FALSE;
  state = FALSE;
  hold = FALSE;
  for (i in start.ind:end.ind) {
    if (hold == TRUE) {
      ret.t = c(ret.t, 1 + stock$RET[i])
    } else {
      ret.t = c(ret.t, 1);
    }
    prev.state = state;
    state = is.good.state(stock, num.clusters, i, look.back);
    hold = state && prev.state;
  }
  return(ret.t);
}

trade.rsr = function(stock, look.back) {
  start.ind = look.back;
  end.ind = length(stock$RET);
  ret.t = c();
  hold = FALSE;
  for (i in start.ind:end.ind) {
    if (hold == TRUE) {
      ret.t = c(ret.t, 1 + stock$RET[i])
    } else {
      ret.t = c(ret.t, 1);
    }
    if (stock$PRC[i] < min(stock$PRC[(i-10):(i-1)])) {
      hold = FALSE;
    } else if (stock$PRC[i] > max(stock$PRC[(i-10):(i-1)])) {
      hold = TRUE;
    }
  }
  return(ret.t);
}

#Organize Data, 
data = read.csv("stats242data_project.csv")
data$DATE = strptime(data$DATE, "%m/%d/%y");
#data$DATE = date.to.ind(data$DATE);
stock.name = c("HAL", "SLB", "EQT", "OXY", 
              "AAPL", "CSCO", "MSFT",
              "BAC", "JPM", "WFC", "C", "STI",
              "ANF", "GES", "SKS", "URBN");
stock.list = get.stock.data.weekly(stock.name, data)

par(mfrow = c(2,2))
#3 Clusters
pick.num.clusters(stock.list[[3]]);
pick.num.clusters(stock.list[[11]]);

#5 Clusters
pick.num.clusters(stock.list[[2]]);
pick.num.clusters(stock.list[[7]]);

#Example stock
num.clusters = 3;
ret = stock.list[[12]]$RET;
ret = ret[2:length(ret)];
init0 = rep(1/num.clusters, num.clusters);
P0 = matrix(1/num.clusters, nrow=num.clusters, ncol=num.clusters);
b0 = get.initial.gaussian(ret, num.clusters);
startval = hmmspec(init=init0, trans=P0, parms.emission=b0, dens.emission=dnorm.hsmm);
h = hmmfit(list(x=ret, N=length(ret)), startval, mstep=mstep.norm, tol=1e-12, maxit=10000);
summary(h)

num.clusters = 3;
look.back = 300;

m.bh = c();
m.hmm = c();
m.rsr = c();
sd.bh = c();
sd.hmm = c();
sd.rsr = c();
sh.bh = c();
sh.hmm = c();
sh.rsr = c();
for (ind in 1:16) {
  #Trade (Buy and Hold)
  ret.bh = trade.buy.and.hold(stock.list[[ind]], look.back);
  port.val.bh = cumprod(ret.bh);
  ret.bh = ret.bh - 1;
  m.bh = c(m.bh, mean(ret.bh));
  sd.bh = c(sd.bh, sd(ret.bh));
  sh.bh = c(sh.bh, mean(ret.bh)/sd(ret.bh));
  
  #Trade (HMM)
  ret.hmm = trade.hmm(stock.list[[ind]], look.back, num.clusters);
  port.val.hmm = cumprod(ret.hmm)
  ret.hmm = ret.hmm - 1;
  m.hmm = c(m.hmm, mean(ret.hmm));
  sd.hmm = c(sd.hmm, sd(ret.hmm));
  sh.hmm = c(sh.hmm, mean(ret.hmm)/sd(ret.hmm));
  
  #Trade (RSR)
  ret.rsr = trade.rsr(stock.list[[ind]], look.back);
  port.val.rsr = cumprod(ret.rsr)
  ret.rsr = ret.rsr - 1;
  m.rsr = c(m.rsr, mean(ret.rsr));
  sd.rsr = c(sd.rsr, sd(ret.rsr));
  sh.rsr = c(sh.rsr, mean(ret.rsr)/sd(ret.rsr));
  
  par(mfrow=c(1,1))
  name=paste0("Portfolio Value Over Time (", stock.name[ind], ")");
  plot(port.val.hmm,col="blue",ylim=c(0, 2.75), type="l", ylab="Portfolio Value", xlab="Time Index", main=name);
  lines(port.val.bh, col="red")
  lines(port.val.rsr, col="green")  
  legend("topleft",  
         c("BH","HMM", "RSR"), 
         lty=c(1,1,1), 
         lwd=c(1,1,1),col=c("red","blue","green"))
  name=paste0("port_val_", stock.name[ind], ".png");
  dev.copy(png,name, width=1000,height=500);
  dev.off();
}
#name=paste0("Price Over Time (", stock.name[[ind]], ")"); 
#plot(stock.list[[1]]$PRC[300:length(stock.list[[ind]]$PRC)], type="l", ylab="Price", xlab="Time Index", main=name)
