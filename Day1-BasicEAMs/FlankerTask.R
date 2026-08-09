drawScreen <- function(txt,txt1=NULL,txt2=NULL,col="black",cex=4,offsetx=0) {
  plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  text(0.5+offsetx,0.5,txt,bty="n",cex=cex,col=col)
  if (!is.null(txt1)) legend("bottomleft",txt1,bty="n",cex=1.25)
  if (!is.null(txt2)) legend("bottomright",txt2,bty="n",cex=1.25)
}


doBlock <- function(type=c("speed","accuracy")[1],trials=sample(rep(1:4,20),80),
                    stim=c("<<<<<",">>>>>",">><>>","<<><<"),stimcol=c("black","black","black","black"),
                    condition=c("congruent","congruent","incongruent","incongruent"),
                    choiceKeys=c("z","m"),responses=c("left","right"),correct=c("left","right","left","right"),
                    offset=c(0,0,0,0),fastrt=0.75,jitter=TRUE) {
  
  
  rrt <- function(choiceKeys=c("m","z")) {
    
    dynamic_readline <- function() {
      while (rstudioapi::isAvailable()) {
        input <- rstudioapi::getConsoleEditorContext()$contents
        
        if (input != "") {
          rstudioapi::sendToConsole("", execute = FALSE)
          return(input)
        }
      }
      readline()
    }
    
    repeat {
      rt <- system.time({
        choice <- dynamic_readline()
      })[3]
      if (choice %in% choiceKeys) break
    }
    setNames(c(choice,rt),c("R","rt"))
  }
  
  dat <- cbind.data.frame(
    E=factor(rep(type,length(trials)),levels=c("accuracy","speed")),
    S=factor(ifelse(trials %in% c(1,3),responses[1],responses[2])),
    # CI=factor(ifelse(trials %in% c(1,2),"congruent","incongruent")),
    CI=factor(condition[trials]),
    R=rep(NA,length(trials)),CR=correct[trials],
    rt=rep(NA,length(trials)),score=numeric(length(trials)))
  
  keyMessage <- paste("Place LEFT index finger on the \"",choiceKeys[1],
                      "\" key and RIGHT index finger on the \"",choiceKeys[2],"\" key")
  keyMessage1 <- paste(choiceKeys[1]," =",responses[1])
  keyMessage2 <- paste(choiceKeys[2]," =",responses[2])
  if (type=="speed") {
    instruct1 <- "100 points for correct"
    instruct2 <- "200 point bonus if also fast"
  } else {
    instruct1 <- "300 points for correct"
    instruct2 <- "100 point penalty for an error"
  }
  drawScreen(paste("Focus on ",type,"\nPress enter to begin"),
             instruct1,instruct2)
  readline()
  drawScreen("Countdown: 3",keyMessage); Sys.sleep(1)
  drawScreen("Countdown: 2",keyMessage); Sys.sleep(1)
  drawScreen("Countdown: 1",keyMessage); Sys.sleep(1)
  drawScreen("",keyMessage); Sys.sleep(1)
  for (i in 1:length(trials)) {
    drawScreen("+",keyMessage1,keyMessage2)
    Sys.sleep(0.5)
    stimi <- stim[trials[i]]
    # add spatial uncertainty
    if (jitter) stimulus <- switch(sample(1:6,1),paste("   ",stimi),paste("  ",stimi),paste("",stimi),
                                   paste(stimi," "),paste(stimi,"  "),paste(stimi,"   ")) else stimulus <- stimi
    drawScreen(stimulus,keyMessage1,keyMessage2,col=stimcol[trials[i]],offsetx=offset[trials[i]])
    dat[i,c("R","rt")] <- rrt(choiceKeys)
    drawScreen(""); Sys.sleep(0.25)
    crct <- responses[choiceKeys==dat$R[i]]==dat$CR[i]
    fast <- as.numeric(dat$rt[i]) < fastrt
    if (type=="speed") {
      if (crct) dat$score[i] <- 100 + 200*fast
      feedback <- paste(ifelse(crct,"Correct: 100 points\n","Error, no points!\n"),
                        ifelse(fast & crct,"SPEED BONUS: 200 points","No speed bonus!"))
      if (crct&fast) col <- "green" else
        if (crct & !fast) col <- "grey" else
          col <-"red"
    } else {
      dat$score[i] <- 300*crct - 100*!crct
      feedback <- paste(ifelse(crct,"Correct: 300 points","Error: -100 points!\n"))
      if (crct) col <- "green" else col <- "red"
    }
    drawScreen(feedback,col=col); Sys.sleep(1)
  }
  drawScreen(paste("Points earned: ",sum(dat$score))); Sys.sleep(4)
  dat$R <- factor(ifelse(dat$R==choiceKeys[1],responses[1],responses[2]))
  dat$rt <- as.numeric(dat$rt)
  dat
}

doExperiment <- function(name,nblocks=10,ntrials=80,choiceKeys=c("z","m"),fastrt=0.7,
                         stim=c("<<<<<",">>>>>",">><>>","<<><<"),stimcol=c("black","black","black","black"),
                         condition=c("congruent","congruent","incongruent","incongruent"),
                         responses=c("left","right"),correct=c("left","right","left","right"),
                         offset=c(0,0,0,0),
                         instruction="Which way does the central arrow point?",jitter=TRUE) {
  par(mfrow=c(1,1))
  if (nblocks %% 2 != 0) stop("nblocks must be even")
  if (ntrials %% 4 != 0) stop("ntrials must be a multiple of 4")
  blocks <- NULL
  drawScreen(paste(instruction,"\n",
                   "Earn as many points as you can!\n",
                   nblocks,"blocks of",ntrials," trials\n Press enter to start"),
             cex=2)
  readline()
  for (i in 1:(nblocks/2))
    blocks <- c(blocks,c("speed","accuracy")[sample(1:2,2)])
  datList <- vector(mode="list",length=nblocks)
  for (i in 1:length(blocks)) {
    datList[[i]] <- doBlock(type=blocks[i],trials=sample(rep(1:4,ntrials/4),ntrials),
                            choiceKeys=choiceKeys,fastrt=fastrt,stim=stim,stimcol=stimcol,condition=condition,
                            responses=responses,correct=correct,jitter=jitter,offset=offset)
  }
  
  data <- cbind.data.frame(subject=factor(rep(name,nblocks*ntrials)),do.call(rbind,datList))
  drawScreen(paste("Experiment done!\nYou earned",sum(data$score)," points!"))
  row.names(data) <- NULL
  if(tolower(stim[1]) %in% colors()){
    data$S_color <- data$S
    data$S <- data$CR
  }
  data
}

