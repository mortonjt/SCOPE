tools = c("SCOPA", "SCOPABW", "CLEAN", "TRIM", "POLY", "TRIMEST", "BASICTOOL")
tools.name = c("Scope", "Scope+bw", "SeqClean", "SeqTrim", "TrimPoly", "TrimEst", "BASIC"); names(tools.name) = tools
colors = c("blue", "black", "red", "yellow", "green", "purple", "orange"); names(colors) = tools
orgs = c("Ara", "Chl", "Human")
orgs.name = c("Arabidopsis", "Chlamydomonus", "Human"); names(orgs.name) = orgs
tech = c("Illumina", "454", "Sanger"); names(tech) = orgs

current_tools = c("SCOPA", "POLY", "TRIMEST", "BASIC")


######### Old 
create.frames <- function(file, obj = "pct_crt", xmin = 0, xmax = Inf, bmin = 0, bmax = Inf, wmin = 0, wmax = Inf, mmin = 0, mmax = Inf) {
  F = read.table(file, header = TRUE);
  F = F[2:length(F[[1]]),];
  
  S1 = F[F$sen >= C$sen & F$x >= xmin & F$x <= xmax & F$b >= bmin & F$b <= bmax & F$w >= wmin & F$w <= wmax & F$m >= mmin & F$m <= mmax, ];
  #S2 = F[F$sen >=  max(F$sen) - 0.02, ];
  S3 = S1[S1[,obj] == max(S1[,obj]), ];
  return(rbind(C, S3))
}

all.frames <- function(obj = "pct_crt", xmin = 0, xmax = Inf, bmin = 0, bmax = Inf, wmin = 0, wmax = Inf, mmin = 0, mmax = Inf) {
  F <- read.table("Ara.T.out", header=TRUE)
  D <- F[1,]
  D$"file" = c('X')
  orgs <- c("Ara", "Chl", "Human")
  types <- c("A", "T")

  for (o in orgs) {
    for (t in types) {
      s <- paste(o, t, "out", sep = ".")
      F <- create.frames(s, xmin=xmin, bmax=bmax)
      F$"file" = s
      D = rbind(D,F)
    }
  }

  return(D[2:length(D[[1]]), c("m","w", "b", "x", "name", "sens", "spec", "pct_crt", "avg_trim", "SoS_trim", "SoS_left", "SoS_right", "file", obj)])
}

check.params <- function(m, w, b, x) {
  F <- read.table("Ara.T.out", header=TRUE)
  D <- F[1,]
  D$file = c('X')
  orgs <- c("Ara", "Chl", "Human")
  types <- c("A", "T")

  for (o in orgs) {
    for (t in types) {
      s <- paste(o, t, "out", sep = ".")
      F <- read.table(s, header = TRUE)
      F$"file" = s
      D = rbind(D, F[1,], F[F$m==m & F$w==w & F$b==b & F$x==x, ])
    }
  }

  return(D[2:length(D[[1]]), c("m","w", "b", "x", "name", "sens", "spec", "pct_crt", "avg_trim", "SoS_left", "SoS_right", "file")])
}



############# Timing 1
plotTiming1 <- function(org, type, usedTools = current_tools, main.msg= FALSE, print.legend= TRUE) {
  F <- read.table(paste(org, type, "timing", "out", sep = "."), header = TRUE)
  F <- F[F$name %in% usedTools,]
  
  plot(c(), c(), xlim = range(F$size), ylim=range(F$run_time), xlab = "Numer of Fragments", ylab="Runtime (s)", main=ifelse(main.msg, main.msg, sprintf("Runtime (%s, %s, poly%s)", org, tech[org], type)))
  for (t in usedTools) {
    X = F[F$name==t, "size"]
    Y = F[F$name==t, "run_time"]
    O = order(X)
    lines(X[O],  Y[O], col = colors[t])
  }

  if (print.legend) {
    legend("topleft", tools.name[usedTools], col = colors[usedTools], lty = c("solid", "solid", "solid"))
  }
}

plotTiming2 <- function(org, type, usedTools = c("SCOPA", "POLY"), main.msg = FALSE) {
  print(paste(org, type, "timing2.1", "out", sep="."))
  F <- read.table(paste(org, type, "timing2.1", "out", sep = "."), header = TRUE)
  
  plot(c(), c(), xlim = range(F$size), ylim=range(F$run_time), xlab = "Numer of Fragments", ylab="Runtime (s)", main=ifelse(main.msg,main.msg,sprintf("Runtime (%s, %s, poly%s)", org, tech[org], type)))
  for (t in usedTools) {
    lines(F[F$name==t, "size"],  F[F$name==t, "run_time"], col = colors[t])
  }

  return(F)
}

  


############## Error robustness
plotError.sen <- function(org, type, usedTools = current_tools, print.legend = TRUE, main.msg=FALSE, y.lim = FALSE) {
  F <- read.table(paste(org, type, "error", "out", sep = "."), header=TRUE)
  F <- F[F$name %in% usedTools,]

  if (length(y.lim) == 1) {
    y.lim = range(F$sens);
  }
  #plot(c(), c(), xlim = range(F$e), ylim=c(0,1), xlab = "Base Call Error Rate", ylab="Sensitivity", main=sprintf("Robustness to Error (Sensitivity)\n(%s, %s, poly%s)", orgs.name[org], tech[org], type))
  plot(c(), c(), xlim = range(F$e), ylim=y.lim, xlab = "Base Call Error Rate", ylab="Sensitivity", main=ifelse(main.msg!=FALSE,main.msg,sprintf("Robustness to Error (Sensitivity)\n(%s, %s, poly%s)", orgs.name[org], tech[org], type)))
  for (t in usedTools) {
    lines(F[F$name==t, "e"],  F[F$name==t, "sens"], col = colors[t])
  }

  if (print.legend) {
    legend("bottomleft", tools.name[usedTools], col = colors[usedTools], lty = c("solid", "solid", "solid"))
  }
}

ScopeVBasic <- function(org, type, field, y.lim = FALSE, p.range = NA) {
    s <- paste(org, type, "svb", "out", sep = ".")
    F <- read.table(s, header = TRUE)

    F <- F[F$name %in% c('BASICTOOL', 'SCOPA'),]

    if (length(y.lim)==1) {
        y.lim = range(F[,field]);
    }

    plot(c(), c(), xlim = range(F$e), ylim = y.lim, xlab = "Base Call Error Rate", main="NEED LABEL")

    
    if (is.na(p)):
      p.range = sort(unique(F$p))
    
    # Print SCOPE
    lines(F[F$name=='SCOPA', 'e'], F[F$name=='SCOPA', field], col = colors['SCOPA'])
    for (p in p.range) {
      print(p)
      lines(F[F$name=='BASICTOOL' & F$p==p,'e'], F[F$name=='BASICTOOL' & F$p==p,field])
    }
}

plotError.correct <- function(org, type, usedTools = current_tools, print.legend = TRUE, main.msg=FALSE, y.lim = FALSE) {
  F <- read.table(paste(org, type, "error", "out", sep = "."), header=TRUE)
  #return(F)
  F <- F[F$name %in% usedTools,]
  
  if (y.lim == FALSE) {
    y.lim = range(F$pct_crt);
  }
  
  plot(c(), c(), xlim = range(F$e), ylim=c(0,1), xlab = "Base Call Error Rate", ylab="% Correct", main=sprintf("Robustness to Error (Correct Trimming)\n(%s, %s, poly%s)", orgs.name[org], tech[org], type))
  #plot(c(), c(), xlim = range(F$e), ylim=y.lim, xlab = "Base Call Error Rate", ylab="% Correct", main=ifelse(main.msg!=FALSE,main.msg,sprintf("Robustness to Error (Correct Trimming)\n(%s, %s, poly%s)", orgs.name[org], tech[org], type)))
  for (t in usedTools) {
    lines(F[F$name==t, "e"],  F[F$name==t, ifelse(type == "A", "pct_crt_left", "pct_crt_right")], col = colors[t])
  }

  if (print.legend) {
    legend("bottomleft", tools.name[usedTools], col = colors[usedTools], lty = c("solid", "solid", "solid"))
  }
}

plotError.both <- function(orgs, types, usedTools = tools, print.legend = 0) {
  quartz(width=10, height = 4*(length(orgs)*length(types)))
  par(mfrow = c(length(orgs)*length(types), 2))

  count = 1
  for (o in orgs) {
    for (t in types) {
      pl = print.legend == count
      print(pl)
        
      plotError.sen(o, t, usedTools, print.legend = pl);
      plotError.correct(o, t, usedTools, print.legend = FALSE);
      count = count + 1
    }
  }
}



plotError.manuscript <- function(usedTools = tools) {
  quartz(width=10, height = 5);
  par(mfrow = c(1,2))

  plotError.sen(c("Chl"), c("A"), print.legend = TRUE, main.msg = "Sensitivity", y.lim=c(0.7, 1))
  plotError.correct(c("Chl"), c("A"), print.legend = FALSE, main.msg = "% correct boundary identification (inner)", y.lim=c(0.4,1));
  #plotError.sen(c("Ara"), c("T"), print.legend = TRUE, main.msg = "Sensitivity, poly(T) tail identification", y.lim=c(0.7,1))
  #plotError.correct(c("Ara"), c("T"), print.legend = FALSE, main.msg = "% correct (inner), poly(T) tail identification", y.lim=c(0.4,1));
}


############################
plot.runTime <- function() {
  quartz(width=5, height=5)
  scale = 10^6
  F <- read.table("Ara.A.timing.2.2.out", header=TRUE)
  B = F$bases/scale
  plot(B, F$run_time, xlab = "Input Size (Mb)", ylab = "Runtime (seconds)", pch = 20, main = "Runtime v. Input Size", axes=FALSE)
  abline(lm(F$run_time ~ B), col = 'red')
  ptsX = pretty(B)
  ptsY = pretty(F$run_time)

  axis(1, at = ptsX, labels=paste(ptsX, "Mb", sep = ""))
  axis(2, at = ptsY, labels=paste(ptsY, "s", sep = ""))

  box()
}
  
plot.trainingSize <- function() {
  quartz(width=5, height=5)
  F <- read.table("output.txt", header=TRUE)


  plot(c(), c(), xlim = range(F$data_size), ylim = range(F$spec), xlab = "Training Set Size", main = "Specificity v. Training Set Size",
       sub = sprintf("Total set size: %d", length(F[1])))
  
  F1 = F[F$retrain == "TRUE", ]
  F2 = F[F$retrain == "FALSE", ]
  lines(F1$data_size, F1$spec, col = 'blue')
  lines(F2$data_size, F2$spec, col = 'green')
}
  
  
################################################################

reload <- function() {
  H00 <<- read.table("Human.A.e0.len.out", header = TRUE)
  H05 <<- read.table("Human.A.e05.len.out", header = TRUE)
  
  H00 <<- H00[H00$tail >= 22 & H00$p.tail <= 270 & H00$tail < 270, c("tail", "p.tail", "name", "sens")]
}

lenPlot <- function(tools = c("SCOPA", "POLY"), F = H00, response = "sens", xlim = NA, ylim = NA, f = function(F) F$p.tail/F$tail, filter = function(F) TRUE, newPlot = TRUE, printLegend = TRUE, xlab = NA, ylab = NA, main = NA, points = TRUE, sub = NA, local.colors = colors) {
  F = F[F$name %in% tools & filter(F), ]
  F$extra = f(F)
  F = F[order(F$extra), ]

  if (newPlot) {
    if (is.na(xlim[1])) xlim = range(F$extra)
    if (is.na(ylim[1])) ylim = range(F[[response]])
    plot(c(), c(), xlim=xlim, ylim=ylim, xlab=xlab, ylab=ifelse(is.na(ylab), response, ylab), main=main, sub=sub)
  }
  
  for (t in tools) {
    G = F[F$name == t, c("extra", response)]
    extraRange = unique(G$extra)
    meanResponse = sapply(extraRange, function(e) mean(G[G$extra == e, response]))
    print(local.colors[t])
    lines(extraRange, meanResponse, xlim=xlim, ylim=ylim, col = local.colors[t])
    #if (points) points(extraRange, meanResponse, xlim=xlim, ylim=ylim, col = colors[t])
  }
  
  if (newPlot & printLegend) legend("right", tools.name[tools], col = colors[tools], lty = "solid")
}

plotRatio1 <- function(F) {
  F = F[F$tail > 0, ]
  M1 = c()
  R = seq(0.25, 10, 0.25)
  
  for (r in R) {
    m = mean(F[F$name == "POLY" & F$p.tail/F$tail == r, "sens"])
    M1 = c(M1, m)
  }
  print(M1)
  
  M2 = c()
  for (r in R) {
    m = mean(F[F$name == "SCOPA" & F$p.tail/F$tail == r, "sens"])
    M2 = c(M2,m)
  }
  print(M2)
  
  plot(c(), c(), xlim = range(R), ylim = range(M1, M2))
  lines(R, M1, col = colors["POLY"])
  lines(R, M2, col = colors["SCOPA"])
}

# p.tail / tail
lenPlot1 <- function(F = DF) {
  lenPlot(F = F, f = function(F) F$p.tail/F$tail, filter = function(F) F$tail > 0 & F$p.tail/F$tail %in% c(0,20,0.5), main = "Ratio: post-tail / tail", xlab = "post-tail / tail", ylab = "sensitivity")
}

# tail / p.tail
lenPlot2 <- function(F = DF) {
  lenPlot(F=F, f = function(F) F$tail/F$p.tail, filter = function(F) F$p.tail > 0 & F$tail/F$p.tail < 10, main = "Ratio: tail / post-tail", xlab = "tail / post-tail", ylab = "sensitivity")
}

# sens. v. tail (averaged over p. tail values)
lenPlot3 <- function(F = DF) {
  lenPlot(F=F, f = function(F) F$tail, main = "tail", xlab = "tail", ylab = "sensitivity")
}

# sens. v. p.tail (avearged over tail values)
lenPlot4 <- function(F = DF) {
  lenPlot(F=F, f = function(F) F$p.tail, main = "p.tail", xlab = "p.tail", ylab = "sensitivity")
}

# sens v. tail for fixed p. tail values
lenPlot5 <- function(ptail.range = seq(10,250,30), xlim = NA, ylim =c(0,1), F=DF, main = NA, sub = NA) {
    lenPlot(F=F, f = function(F) F$tail, filter = function(F) F$p.tail == ptail.range[1], xlim = xlim, ylim = ylim, main = ifelse(is.na(main), "variable t with fixed p", main), xlab = "tail length",ylab = "sensitivity", sub = sub)
    for (i in ptail.range[-1]) {
      lenPlot(F=F, f = function(F) F$tail, filter = function(F) F$p.tail == i, newPlot = FALSE)
    }
}

# sens v. p. tail for fixed tail values
lenPlot6 <- function(tail.range = seq(30,260,10), xlim = NA, ylim = c(0,1), F=DF, main = NA) {
  lenPlot(F=F, f = function(F) F$p.tail, filter = function(F) F$tail == tail.range[1], xlim = xlim, ylim = ylim, main = ifelse(is.na(main), "variable p with fixed t", main), xlab = "p.tail",ylab = "sensitivity")
  for (i in tail.range[-1]) {
    lenPlot(F=F, f = function(F) F$p.tail, filter = function(F) F$tail == i, newPlot = FALSE)
  }
}

tail.v.ptail <- function(F = DF, tools = c("SCOPA", "POLY"), sens.min = 0.98, sens.max = 1) {
  F = F[F$sens >= sens.min & F$sens <= sens.max & F$name %in% tools, ]
  
  plot(c(), c(), xlim = range(F$tail), ylim = range(F$p.tail), xlab = "tail", ylab = "ptail", main = sprintf("ptail v tail at fixed sensitivity range: %5.2f-%5.2f", sens.min, sens.max))
  for (t in tools) {
    points(F[ F$name == t, c("tail", "p.tail") ], col = colors[t])
  }
}
  

#######################################
# Set of specific plots for analysis of the (1-270)x(1-270) data set.


#####################################
# For manuscript
 
# Plot Sensitivity agaist post-tail length / tail length
Figure1 <- function(F = H00) {
  F = F[F$tail > 0, ]
  M1 = c()
  R = seq(0.25, 10, 0.25)
  
  for (r in R) {
    m = mean(F[F$name == "POLY" & F$p.tail/F$tail == r, "sens"])
    M1 = c(M1, m)
  }
  
  M2 = c()
  for (r in R) {
    m = mean(F[F$name == "SCOPA" & F$p.tail/F$tail == r, "sens"])
    M2 = c(M2,m)
  }
  
  plot(c(), c(), xlim = range(R), ylim = range(M1, M2), ylab = "Sensitivity", xlab = "post-tail length / tail length", main = "Sensitivity v. (post-tail length) / (tail length) Ratio", sub="(a)")
  lines(R, M1, col = colors["POLY"])
  lines(R, M2, col = colors["SCOPA"])
  usedTools = c("SCOPA", "POLY")
  legend("bottomleft", tools.name[usedTools], col = colors[usedTools], lty = c("solid", "solid", "solid"))
}

Figure2 <- function(F = H00) {
  #lenPlot5(seq(50, 250, 50), main = "Sensitivity v. tail length at fixed values of post-tail length", sub = "(b)")
  ptail.range = seq(50,250,50)
  F = F[F$p.tail %in% ptail.range,]
  xrange = range(F$tail)
  yrange = range(F$sens)
  plot(c(), c(), xlim = xrange, ylim = yrange, xlab = "Tail length", ylab = "sensitivity", main = "Sensitivity v. Tail Length at fixed values of\npost-tail segment length", sub="(b)")

  tmp.colors = c("green", "red", "yellow", "purple", "black")
  # First: plot scope

  for (i in 1:length(ptail.range)) {
    C = c(tmp.colors[i])
    names(C) = c("POLY")
    lenPlot(tools = c("POLY"), F = F, response = "sens", f = function(F) F$tail, filter = function(F) F$p.tail == ptail.range[i], newPlot = FALSE, printLegend = FALSE, local.colors = C)
  }
  lenPlot(tools = c("SCOPA"), F = F, response = "sens", f = function(F) F$tail, filter = function(F) TRUE, printLegend = FALSE, newPlot = FALSE)
  legend("bottom", c("SCOPA (all post-tail lenths)", "POLY: post-tail length = 50",  "POLY: post-tail length = 100",  "POLY: post-tail length = 150",  "POLY: post-tail length = 200", "POLY: post-tail length = 250"), col = c("blue", "green", "red", "yellow", "purple", "black"), lty = replicate(5, "solid"), bg = "white")
}
  

Figure <- function(F = H00, file = NA) {
  height = 5
  width = 12
  if (is.na(file)) {
    quartz(height=height, width=width)
  }
  else {
    pdf(file, width=width, height=height)
  }
  par(mfrow = c(1,2))
  Figure1()
  Figure2()
  if (!is.na(file)) {
    dev.off()
  }
}
