arrayMissPattern<-function()
{
   require(RGtk2)
   require(cairoDevice)
   require(spatstat)
   
   Result.Save<<-NULL
   high.list<<-NULL
   high.chip.list<<-NULL
   gene.location<<-NULL
   true.flag<<-FALSE

#############################################################################################
## control GeneList, ChipList, and output Table

geneListSelectRow.handler <- function(w, r, c, e, u = NULL)
{     high.list<<-sort(cbind(high.list,(r+1)))
      high.list<<-unique(high.list)
}

geneListUnselectRow.handler <- function(w, r, c, e, u = NULL)
{     high.list<<-high.list[high.list!=(r+1)]
      high.list<<-unique(high.list)
}

chipListSelectRow.handler <- function(w, r, c, e, u = NULL)
{    high.chip.list<<-sort(cbind(high.chip.list,(r+1)))
      high.chip.list<<-unique(high.chip.list)
}

chipListUnselectRow.handler <- function(w, r, c, e, u = NULL)
{     high.chip.list<<-high.chip.list[high.chip.list!=(r+1)]
      high.chip.list<<-unique(high.chip.list)
}

showMissingList<- function(missing.Name)
{    if(ncol(missing.Name)==3)
     {  temp.p<-9
        temp<-matrix(rep(" ",temp.p*nrow(missing.Name)),ncol=temp.p)
        missing.Name<-cbind(missing.Name,temp)
     }
     for(i in 1: nrow(missing.Name))
        gtkCListAppend(geneList,as.character(missing.Name[i,]))
}

showChipList<- function()
{   chip.name<<-as.matrix(colnames(data.missing))
     for(i in 1: nrow(chip.name))
        gtkCListAppend(ChipList,as.character(chip.name[i,]))
}

removeMissingList <- function(n)
{   for(i in (n-1):0)
       gtkCListRemove(geneList,i)
}

warning.gtk <- function(messages) 
{  warn.win <- gtkWindow(show = FALSE)
    box <- gtkVBox(TRUE, 0)
    box$SetUsize(-1, 80)
    box$PackStart(gtkLabel(messages))
    box$PackStart(box.small <- gtkHBox(FALSE, 0), expand = TRUE, fill = FALSE)
    box.small$PackStart(button <- gtkButton("Ok"), expand = TRUE, fill = FALSE)
    button$SetUsize(60, 25)
    button$AddCallback("clicked", function(x) {
      warn.win$Hide()
      warn.win$Destroy()
    })
    warn.win$Add(box)
    warn.win$Show()
}
  
############################################################################################# 
### In File menu
## Read microaray data with missing

Mnu.read.file.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open Microarray data file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         data.missing<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
         rownames(data.missing)<<-as.character(data.missing[,1])
         data.missing<<-data.missing[,-1]
         missing.indicator(data.missing)
         missing.pattern<<-matrix(as.numeric(is.na(data.missing)),nrow=nrow(data.missing))

         showChipList()
         contents <- paste(" Read ", exp.info.name, sep=" ")
         contents <- paste(contents,"\n",sep=" ")
         gtkTextBufferInsertAtCursor(new.text,contents)
         gtkTextViewSetBuffer(outputText,new.text)

         win$Hide()
         win$Destroy()
    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })
}

## Read Gene Location with pin information

Mnu.read.grid.file.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open Gene Grid Information file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         gene.location<<-NULL
         gene.location<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
         x<<-gene.location[,1]
         y<<-gene.location[,2]
         contents <- paste(" Read ", exp.info.name, sep=" ")
         contents <- paste(contents,"\n",sep=" ")
         gtkTextBufferInsertAtCursor(new.text,contents)
         gtkTextViewSetBuffer(outputText,new.text)
         win$Hide()
         win$Destroy()

    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })
}

Mnu.read.block.file.handler<-function(w,u=NULL)
{  win<- gtkFileSelection("Open Gene Block Information file")
   ok <- win[["OkButton"]]
   cancel <- win[["CancelButton"]]
   ok$AddCallback("clicked", function(x) {
         gene.block<<-NULL
         gene.block<<-read.csv(exp.info.name<<-win$GetFilename(),header=T)
         contents <- paste(" Read ", exp.info.name, sep=" ")
         contents <- paste(contents,"\n",sep=" ")
         gtkTextBufferInsertAtCursor(new.text,contents)
         gtkTextViewSetBuffer(outputText,new.text)
         win$Hide()
         win$Destroy()

    }) 
   cancel$AddCallback("clicked", function(x) {
         win$Hide()
         win$Destroy()
   })
}


## Exit

Mnu.exit.handler<-function(w,u=NULL)
{  main$Hide()
   main$Destroy()
   dev.off()

}

missing.indicator<-function(data.missing)
{  missing.I<<-NULL
   missing.Name<<-NULL
   gene.name<-rownames(data.missing)
   chip.name<-colnames(data.missing)
   for(k in 1:ncol(data.missing))
   {  temp.I<-which(is.na(data.missing[,k]))
      missing.Name<<-rbind(missing.Name,cbind(rep(chip.name[k],length(temp.I)),gene.name[temp.I]))
      missing.I<<-rbind(missing.I,cbind(rep(k,length(temp.I)),temp.I))
   }
   colnames(missing.Name)<<-c("Chip","Gene")
   colnames(missing.I)<<-c("Chip","Gene")
   missing.Name.temp<<-missing.Name
   miss.temp<<-missing.I
}

#############################################################################################
## in Exploring missing patterns menu
## missing rates per chip
Mnu.Missing.rates.chip.handler<-function(w,u=NULL)
{  Cairo()
   a<-apply(data.missing,2,function(x){sum(is.na(x))})/nrow(data.missing)
   barplot(a,ylim=c(0,min(1,2*max(a))),main="missing rates for each chip",las=2,cex.axis=par(0.7))  
}

## missing rates per gene
Mnu.Missing.rates.gene.handler<-function(w,u=NULL)
{  Cairo()
   a<-apply(data.missing,1,function(x){sum(is.na(x))})/ncol(data.missing)
   plot(1:length(a),a,type='h',main="missing rates of each gene",xlab="gene",ylab="missing rate")
}

## missing rates per block
Mnu.Missing.rates.block.handler<-function(w,u=NULL)
{  if(is.null(gene.block))
   {  warning.gtk("Need PIN(block) information of genechip")
   } else
   {  missing.pattern<-matrix(as.numeric(is.na(data.missing)),nrow=nrow(data.missing))
      if(ncol(gene.block)==2)
      { pin<-paste(gene.block[,1],gene.block[,2],sep="")
      } else
      { pin<-gene.block
      }
      Cairo()
      a<-apply(missing.pattern,1,sum)/ncol(data.missing)
      boxplot(a ~factor(pin),las=2,main="boxplot for missing rates of each gene",xlab="block")
   }
}

## Overall missing pattern
Mnu.Missing.pattern1.handler<-function(w,u=NULL)
{  Cairo()
   a<-data.missing
   a[!is.na(data.missing)]<-10
   a[is.na(data.missing)]<-0
   image(t(as.matrix(a)),axes=F,ylab="genes")
   axis(1,at=seq(from=0,to=1,length.out=ncol(data.missing)),lab=colnames(data.missing),lty=2,las=2)
}

## Missing pattern for location in gene chip
Mnu.Missing.pattern2.handler<-function(w,u=NULL)
{  if(is.null(gene.block))
   {  warning.gtk("Need location information of genechip")
   } else
   {  Cairo()
       missing.count<-apply(missing.pattern,1,sum) 
       n<-length(table(x))
       p<-length(table(y))
       a.image<-matrix(rep(1,n*p),ncol=p)
       missing.count[missing.count>round(ncol(data.missing)/2)]<-round(ncol(data.missing)/2)
       for(i in 1:length(missing.count))
           a.image[gene.location[i,1],gene.location[i,2]]<-missing.count[i]
       image(as.matrix(a.image),axes=F,ylab="y.location",xlab="x.location")
    }
}

## Exploring Missing Pattern for each chip

Mnu.Missing.pattern.chip.handler<-function(w,u=NULL)
{  if(length(high.chip.list)!=1)
   {   warning.gtk("Please select one chip at a time")
   } else  if(is.null(gene.block))
   {  warning.gtk("Need location information of genechip")
   } else
   {   if(is.null(gene.location))
       {  warning.gtk("Need locatio  information of genechip")
       } else
       {  Cairo()
         
          missing.count<-as.numeric(is.na(data.missing[,high.chip.list]))
          n<-length(table(x))
          p<-length(table(y))
          a.image<-matrix(rep(1,n*p),ncol=p)
          for(i in 1:length(missing.count))
             a.image[x[i],y[i]]<-1-missing.count[i]

          image(as.matrix(a.image),axes=F,ylab="y.location",xlab="x.location",
               main=as.character(chip.name[high.chip.list,1]))
      }	        
  }
}

## Intensity based exploration

Mnu.Missing.pattern.intensity.handler<-function(w,u=NULL)
{  Cairo()
   mean.intensity<-apply(data.missing,1,mean,na.rm=T)
   missing.pattern<-matrix(as.numeric(is.na(data.missing)),nrow=nrow(data.missing))
   a<-apply(missing.pattern,1,sum)/ncol(data.missing)
   plot(mean.intensity,a,ylab="missing.rate")
}

Mnu.Ripley.K.handler<-function(w,u=NULL)
{   print("Ripley.K")
    missing.n.chip<-apply(data.missing,2,function(x){sum(as.numeric(is.na(x)))})
    missing.pattern<-matrix(as.numeric(is.na(data.missing)),nrow=nrow(data.missing))

    Lm<-NULL; Ls<-NULL
    T<-data.missing[,1]
    x<-NULL;  y<-NULL
    for(i in 1:nrow(missing.pattern))
    {   if(is.na(T[i]))
        {   x<-c(x,gene.location[i,1]/max(gene.location[,1]))
            y<-c(y,gene.location[i,2]/max(gene.location[,2]))
        }
    }
    dd<-list(x=x,y=y)
    X<-as.ppp(dd,c(0,1,0,1))
    K <- Kest(X, correction=c("isotropic"))

# isotropic
    Cairo()
 
    plot(K$r,sqrt(K$theo/pi)-K$r,type='l',xlim=c(0,0.25),ylim=c(-0.1,0.25),xlab="t",ylab="L(t)-t",col="blue",lwd=3)
    Lm.temp<-max(abs(sqrt(K$iso/pi)-K$r))
    if(Lm.temp>1.42/sqrt(missing.n.chip[1]))
    { lines(K$r,sqrt(K$iso/pi)-K$r,lty=1,col="red")
    } else
    {lines(K$r,sqrt(K$iso/pi)-K$r,lty=1,col="grey")}

    Lm<-rbind(Lm,c(1,max(abs(sqrt(K$iso/pi)-K$r))))
    Ls<-rbind(Ls,c(1,sum(abs(sqrt(K$iso/pi)-K$r))))

    for(ii in 2:ncol(data.missing))
   {  T<-data.missing[,ii]
       x<-NULL
       y<-NULL
       for(i in 1:nrow(missing.pattern))
       {    if(is.na(T[i]))
            {   x<-c(x,gene.location[i,1]/max(gene.location[,1]))
                 y<-c(y,gene.location[i,2]/max(gene.location[,2]))
            }
       }
       dd<-list(x=x,y=y)
       X<-as.ppp(dd,c(0,1,0,1))
       K <- Kest(X, correction=c("border","isotropic"))
       Lm.temp<-max(abs(sqrt(K$iso/pi)-K$r))
       if(Lm.temp>1.42/sqrt(missing.n.chip[ii]))
       { lines(K$r,sqrt(K$iso/pi)-K$r,lty=1,col="red")
       } else
       {lines(K$r,sqrt(K$iso/pi)-K$r,lty=1,col="grey")}
        Lm<-rbind(Lm,c(ii,max(abs(sqrt(K$iso/pi)-K$r))))
        Ls<-rbind(Ls,c(ii,sum(abs(sqrt(K$iso/pi)-K$r))))
    }

     CV.5<-1.42/sqrt(missing.n.chip[Lm[,1]])
     Cairo()
     plot(Lm,type='n',xlab="chip",ylab="Lm",main="")
     title(main="Lm=sup(|L(t)-t|)")
     for(i in 1:nrow(Lm))
     {   if(Lm[i,2]>1.42/sqrt(missing.n.chip[i]))
          {   lines(c(Lm[i,1],Lm[i,1]),c(0,Lm[i,2]),col="red")
          } else
          { lines(c(Lm[i,1],Lm[i,1]),c(0,Lm[i,2]),col="grey")}
     }
     L.result<-cbind(Lm,CV.5)  
     contents <- paste(" Chip   ", "  Lm", "    CV(0.05)", "\n",sep=" ")
     for(i in 1:nrow(L.result))
     {    contents <- paste(contents, chip.name[i,1],round(L.result[i,2],4),round(L.result[i,3],4),"\n",sep="    ")
     }
      gtkTextBufferInsertAtCursor(new.text,contents)
      gtkTextViewSetBuffer(outputText,new.text)

}

#############################################

############################################
############################################
# Main GUI
############################################
###########################################


 #-----  Main Window  -----#
   main <<- gtkWindow (show = FALSE) 
   main$SetTitle("arrayMissPattern")
   gtkWindowSetDefaultSize(main,460,300)
#   gtkWindowSetPolicy(main,TRUE,TRUE,FALSE)  #allow top level window resizable
   gtkWindowSetResizable(main,TRUE)
   main$SetUposition(40,100) #gtkWidgetSetUposition



 #-----  Chip List  -----#     (Temporary)
   ChipList <<- gtkCListNew(1)
#   gtkCListSetColumnWidth(ChipList,0,85)

   gtkCListColumnTitlesShow(ChipList)
   label11 <- gtkLabelNew("Chip name")	
   gtkCListSetColumnWidget(ChipList,0,label11)
   gtkCListColumnTitlesPassive(ChipList)
   gtkCListSetColumnResizeable(ChipList,0,FALSE)
   gtkCListSetSelectionMode(ChipList,GtkSelectionMode[4])

 #----- output screen -----#
    outputText <<- gtkTextViewNew()
    gtkTextViewSetEditable(outputText,FALSE)
#    outputText$SetEditable(FALSE)
    new.text<<-gtkTextBufferNew()  

#----- Scroll Window 1&2 -----#
 
 
   scrollwindow2 <- gtkScrolledWindowNew()	# to go into the geneList frame
   scrollwindow2$SetUsize(300,350)
   gtkScrolledWindowSetPolicy(scrollwindow2,GtkPolicyType[2],GtkPolicyType[2])
   scrollwindow2$Add(outputText)
   gtkContainerSetBorderWidth(scrollwindow2,3)

   scrollwindow3 <- gtkScrolledWindowNew()	# to go into the geneList frame
   scrollwindow3$SetUsize(150,100)
   gtkScrolledWindowSetPolicy(scrollwindow3,GtkPolicyType[2],GtkPolicyType[2])
   scrollwindow3$Add(ChipList)
   gtkContainerSetBorderWidth(scrollwindow3,3)

#-----------------------------------#
#               menu
#-----------------------------------#

   menu1 <- gtkMenuNew()
   menu2 <- gtkMenuNew()
   menu3 <- gtkMenuNew()
   menu4 <- gtkMenuNew()
   menu21 <- gtkMenuNew()
   menu41 <- gtkMenuNew()
   menu42 <- gtkMenuNew()
 
#----------

   show1 <- gtkMenuItemNewWithLabel("File")
   item11 <- gtkMenuItemNewWithLabel("Read Microarray data with missing")
   item12 <- gtkMenuItemNewWithLabel("Read GeneChip Location information")
   item13 <- gtkMenuItemNewWithLabel("Read GeneChip pin information")
   item14 <- gtkMenuItemNewWithLabel("Exit")

# buttonNames<-c("kNNimpute"","BPCAimpute","LSimpute","LLSimpute","RLSPimpute","SVRimpute")
   show2 <- gtkMenuItemNewWithLabel("Explore Missing Patterns")   
   item21 <- gtkMenuItemNewWithLabel("Missing Rates")
   item211 <- gtkMenuItemNewWithLabel("per chip")
   item212 <- gtkMenuItemNewWithLabel("per gene")
   item213 <- gtkMenuItemNewWithLabel("per pin")
   item22 <- gtkMenuItemNewWithLabel("Overall Missing Pattern")
   item23 <- gtkMenuItemNewWithLabel("Missing Pattern for location of chip")
   item24 <- gtkMenuItemNewWithLabel("Explore Missing Pattern for each chip")
   item25 <- gtkMenuItemNewWithLabel("Intensity Based Exploration")

   show3 <- gtkMenuItemNewWithLabel("Test for Random Pattern")   
   item31 <- gtkMenuItemNewWithLabel("Ripley's K")
 
  
#-----------

   gtkMenuShellAppend(menu1,item11)
   gtkMenuShellAppend(menu1,item12)
   gtkMenuShellAppend(menu1,item13)
   gtkMenuShellAppend(menu1,item14)
 
   gtkMenuShellAppend(menu2,item21) # missing rate
   gtkMenuShellAppend(menu2,item22) # missing pattern1
   gtkMenuShellAppend(menu2,item23) # missing pattern 2 
   gtkMenuShellAppend(menu2,item24) # missing pattern 2 
   gtkMenuShellAppend(menu2,item25) # missing pattern 2 

   gtkMenuShellAppend(menu21,item211) # per chip
   gtkMenuShellAppend(menu21,item212) # per gene
   gtkMenuShellAppend(menu21,item213) # per pin

   gtkMenuShellAppend(menu3,item31) # setting parameters
  

#------------

   gtkWidgetShow(show1)
   gtkWidgetShow(item11) 
   gtkWidgetShow(item12) 
   gtkWidgetShow(item13) 
   gtkWidgetShow(item14) 
 
   gtkWidgetShow(show2)    
   gtkWidgetShow(item21) 
   gtkWidgetShow(item22) 
   gtkWidgetShow(item23) 
   gtkWidgetShow(item24) 
   gtkWidgetShow(item25) 
   gtkWidgetShow(item211) 
   gtkWidgetShow(item212) 
   gtkWidgetShow(item213) 

   gtkWidgetShow(show3)    
   gtkWidgetShow(item31) 
   
#-----------

   gtkMenuItemSetSubmenu(show1,menu1)
   gtkMenuItemSetSubmenu(show2,menu2)
   gtkMenuItemSetSubmenu(item21,menu21)
   gtkMenuItemSetSubmenu(show3,menu3)  
# --------

   menuBar<-gtkMenuBarNew()
   gtkWidgetShow(menuBar)

   gtkMenuShellAppend(menuBar,show1)
   gtkMenuShellAppend(menuBar,show2)
   gtkMenuShellAppend(menuBar,show3)


#------ callback function

   gtkAddCallback(item11, "activate", Mnu.read.file.handler)
   gtkAddCallback(item12, "activate", Mnu.read.grid.file.handler)
   gtkAddCallback(item13, "activate", Mnu.read.block.file.handler)
   gtkAddCallback(item14, "activate", Mnu.exit.handler)

   gtkAddCallback(item211, "activate", Mnu.Missing.rates.chip.handler)
   gtkAddCallback(item212, "activate", Mnu.Missing.rates.gene.handler)
   gtkAddCallback(item213, "activate", Mnu.Missing.rates.block.handler)
   gtkAddCallback(item22, "activate", Mnu.Missing.pattern1.handler)
   gtkAddCallback(item23, "activate", Mnu.Missing.pattern2.handler)
   gtkAddCallback(item24, "activate", Mnu.Missing.pattern.chip.handler)
   gtkAddCallback(item25, "activate", Mnu.Missing.pattern.intensity.handler)
 
   gtkAddCallback(item31, "activate", Mnu.Ripley.K.handler)
   
   gtkAddCallback(ChipList, "select-row", chipListSelectRow.handler)
   gtkAddCallback(ChipList, "unselect-row", chipListUnselectRow.handler)

#------------------------------------------------------------------#
#------------------------------------------------------------------#


#-----  Frames  -----#  
#-----  Frame : GeneList List -----#
 

#----- Frame : OutputTable ---------#
   outputTableFrame<-gtkFrameNew("Results")
   gtkContainerAdd(outputTableFrame,scrollwindow2)
   gtkWidgetShow(outputTableFrame)


#----- Frame : ChipList ---------#
   ChipListFrame<-gtkFrameNew("Chip List")
   gtkContainerAdd(ChipListFrame,scrollwindow3)
   gtkWidgetShow(ChipListFrame)


#----- Packing and Show -----#

    Hpan <-gtkHPanedNew(FALSE) # to put GeneListFrame & ChipListFrame together
    Vpan <-gtkVPanedNew(FALSE)
   gtkPanedPack1(Hpan,ChipListFrame,TRUE,TRUE) 
   gtkPanedPack2(Hpan,outputTableFrame,TRUE,TRUE) 
   Hpan$Show() 
   mainLayout <- gtkVBoxNew(FALSE,0)
   gtkBoxPackStart(mainLayout,menuBar,FALSE,FALSE,0)
   gtkBoxPackStart(mainLayout,Hpan,TRUE,TRUE,0)
   main$Add(mainLayout)
   main$Show()

}
