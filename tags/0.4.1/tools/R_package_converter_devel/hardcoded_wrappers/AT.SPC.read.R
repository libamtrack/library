# R function for reading of spc files 
#
# started 2010/02/18, sgre
# as R script
#
# revised 2010/05/18, fk
# added: show.bin/show.bin.width feature, set type of mean.
# Bug fixed: tmp.data$Cum changed to tmp.data$HCum after tmp.data$H.
#
# revised 2010/07/28, sgre
# removed DDD handling
#
# revised 2010/08/04, sgre
# now choice of endian
#
# revised 2010/11/02, sgre
# added to R package


AT.SPC.read <- function( file.name, 
                         endian         = c("big", "little")[1], 
                         mean           = c(geometric = 0, arithmetic = 1)[2],
                         raw            = FALSE,
                         compress       = TRUE)
{
    ###########################
    # DEFINE INTERNAL FUNCTIONS
    ###########################

    # 1. function for reading tags
    read.tag    <-    function(    file.handle, endian){
        tag         <-    readBin(file.handle, integer(), endian = endian)
        length      <-    readBin(file.handle, integer(), endian = endian)
        return(list(    tag       =    tag,
                        length    =    length))
    }

    # 2. tag codes
    SPC.tags    <-    data.frame(    tag        =    1:20,    
                        TRiP.name    =    c(    "FILETYPE",            "FILEVERSION",            "FILEDATE",            "TARGNAME",                    "PROJNAME",
                                        "B",                "P",                    "N",                "NZ",                        "Z",
                                        "N",                "NS",                    "S",                "CUM",                    "NC",
                                        "NE",                "E",                    "EREF",            "HISTO",                    "RUNNINGSUM"),
                        real.name    =    c(    "file.type",        "file.version",            "file.date",        "target.name",                "projectile.name",
                                        "beam.energy.MeV.u",    "peak.position.g.cm2",        "normalization",        "number.of.depth.steps",        "depth.g.cm2",
                                        "normalization",        "number.of.particle.species",    "Z.and.A",            "cumulated.number.of.fragments",    "not.used",
                                        "energy.MeV.u",        "number.of.energy.bins",    "EREF",            "number.of.particles.per.bin",    "cum.number"),
                        format    =    c(    "String",            "String",                "String",            "String",                    "String",
                                        "double",            "double",                "double",            "ulong",                    "double",
                                        "double",            "ulong",                "2double2long",        "double",                    "ulong",
                                        "ulong",            "double",                "ulong",            "double",                    "double"))


    SPC.format.lookup    <-    function(SPC.table, tag){
        return(as.character(SPC.table$format[match(tag,SPC.table$tag)]))
    }


    SPC.TRiP.name.lookup    <-    function(SPC.table, tag){
        return(as.character(SPC.table$TRiP.name[match(tag,SPC.table$tag)]))
    }

    SPC.real.name.lookup    <-    function(SPC.table, tag){
        return(as.character(SPC.table$real.name[match(tag,SPC.table$tag)]))
    }

    # 6. function for reading data items
    read.item    <-    function(    file.handle, endian){
        # read tag first
        item              <-    read.tag(    file.handle, endian)
        item$format       <-    SPC.format.lookup(    SPC.tags, item$tag)
        item$TRiP.name    <-    SPC.TRiP.name.lookup(    SPC.tags, item$tag)
        item$real.name    <-    SPC.real.name.lookup(    SPC.tags, item$tag)

        # branch according to tag
        if(item$format == "String"){
            item$data <- rawToChar( readBin( file.handle, raw(), n = item$length, endian = endian))    
        }

        if(item$format == "double"){
            item$data <- readBin( file.handle, double(), size = 8, n = floor(item$length/8), endian = endian)    
        }

        if(item$format == "2double2long"){
            Z         <- readBin( file.handle, double(),  size = 8, n = 1, endian = endian)
            A         <- readBin( file.handle, double(),  size = 8, n = 1, endian = endian)
            lZ        <- readBin( file.handle, integer(), size = 4, signed = TRUE, n = 1, endian = endian)
            lA        <- readBin( file.handle, integer(), size = 4, signed = TRUE, n = 1, endian = endian)
            item$data <- list( Z  = Z,
                               A  = A,
                               lZ = lZ,
                               lA = lA)
        }

        if(item$format == "ulong"){
            item$data    <-    readBin(    file.handle, integer(), size = 8, signed = FALSE, n = floor(item$length/8), endian = endian)    
        }
        return(item)
    }

    #######################
    # START ACTUAL FUNCTION
    #######################

    # Open file for binary access
    to.read                 <- file(file.name, "rb")
    # read file header
    spc                     <- list()
    spc$file.type           <- read.item( to.read, endian)$data
    spc$file.version        <- read.item( to.read, endian)$data
    spc$file.date           <- read.item( to.read, endian)$data
    spc$target.name         <- read.item( to.read, endian)$data
    spc$projectile.name     <- read.item( to.read, endian)$data
    spc$beam.energy.MeV.u   <- read.item( to.read, endian)$data
    spc$peak.position.g.cm2 <- read.item( to.read, endian)$data
    spc$normalization       <- read.item( to.read, endian)$data
    spc$n.depth.steps       <- read.item( to.read, endian)$data

    spc$data                <- NULL

    for (i in 1:spc$n.depth.steps){
        #i<-1
        depth.g.cm2             <-    read.item( to.read, endian)$data
        normalization           <-    read.item( to.read, endian)$data
        n.particle.species      <-    read.item( to.read, endian)$data
        
        for (j in 1:n.particle.species){
            #j<-1
            #j<-3
            Z.and.A                 <-    read.item( to.read, endian)$data
            Cum                     <-    read.item( to.read, endian)$data
            NC                      <-    read.item( to.read, endian)$data
            n.energy.bins           <-    read.item( to.read, endian)$data
            
            tmp.data                <-    data.frame(   depth.step          =    rep(i,                     n.energy.bins),
                                                        depth.g.cm2         =    rep(depth.g.cm2,           n.energy.bins),
                                                        n.particle.species  =    rep(n.particle.species,    n.energy.bins),
                                                        particle.species    =    rep(j,                     n.energy.bins),
                                                        Z                   =    rep(Z.and.A$Z,             n.energy.bins),
                                                        A                   =    rep(Z.and.A$A,             n.energy.bins),
                                                        lZ                  =    rep(Z.and.A$lZ,            n.energy.bins),
                                                        lA                  =    rep(Z.and.A$lA,            n.energy.bins),
                                                        Cum                 =    rep(Cum,                   n.energy.bins),
                                                        n.energy.bins       =    rep(n.energy.bins,         n.energy.bins),
                                                        E.low.MeV.u         =    rep(0.0,                   n.energy.bins),
                                                        E.mid.MeV.u         =    rep(0.0,                   n.energy.bins),
                                                        E.high.MeV.u        =    rep(0.0,                   n.energy.bins),
                                                        H                   =    rep(0.0,                   n.energy.bins),
                                                        HCum                =    rep(0.0,                   n.energy.bins))
                                        
            tmp.item                <-    read.item( to.read, endian)

            if(tmp.item$tag == 18){    #EREFCOPY        
                ref.species             <-    tmp.item$data
                ii                      <-    spc$data$depth.step == i & spc$data$particle.species == (ref.species + 1)
                tmp.data$E.low.MeV.u    <-    spc$data$E.low.MeV.u[ii]
                tmp.data$E.high.MeV.u   <-    spc$data$E.high.MeV.u[ii]
                tmp.data$E.mid.MeV.u    <-    spc$data$E.mid.MeV.u[ii]
            }else{
                E.bins.MeV.u            <-    tmp.item$data
                tmp.data$E.low.MeV.u    <-    E.bins.MeV.u[-length(E.bins.MeV.u)]
                tmp.data$E.high.MeV.u   <-    E.bins.MeV.u[-1]
                # calculate E.mid.MeV.u 
                if(mean == 0){
                    tmp.data$E.mid.MeV.u    <-    sqrt(tmp.data$E.low.MeV.u * tmp.data$E.high.MeV.u)
                }
                if(mean == 1){
                    tmp.data$E.mid.MeV.u    <-    (tmp.data$E.low.MeV.u + tmp.data$E.high.MeV.u)/2
                }            
            }

            tmp.data$H              <-    read.item( to.read, endian)$data
            tmp.data$HCum           <-    read.item( to.read, endian)$data[-1]
            
            if(is.null(spc$data)){
                spc$data                <-    tmp.data
            }else{
                spc$data                <-    rbind.data.frame(spc$data,tmp.data)
            }
            #cat(paste("Read particle species", j, "of", n.particle.species, "(Z = ", Z.and.A$lZ, "A = ", Z.and.A$lA, ")\n"))
        }
        cat(paste("Read depth step", i, "of", spc$n.depth.steps, "\n"))
    }

    # close file
    close(to.read)

    # add bin width
    spc$data$DE.MeV.u    <- spc$data$E.high.MeV.u - spc$data$E.low.MeV.u

    # Remove zero fluence bins if requested
    if(compress){
      ii         <-  spc$data$H != 0
      spc$data   <-  spc$data[ii,]
    }

    # reformat if not requested otherwise
    if(!raw){
       spc$data$particle.no <- AT.particle.no.from.Z.and.A( Z = spc$data$Z,
                                                            A = spc$data$A)$particle.no
       spc$data$Z                  <- NULL
       spc$data$A                  <- NULL
       spc$data$lZ                 <- NULL
       spc$data$lA                 <- NULL
       spc$data$n.particle.species <- NULL
       spc$data$particle.species   <- NULL
       spc$data$Cum                <- NULL
       spc$data$n.energy.bins      <- NULL
       spc$data$E.low.MeV.u        <- NULL
       spc$data$E.high.MeV.u       <- NULL
       spc$data$HCum               <- NULL

       spc$data$E.MeV.u            <- spc$data$E.mid.MeV.u
       spc$data$E.mid.MeV.u        <- NULL

       spc$data$fluence.cm2        <- spc$data$H * spc$data$DE.MeV.u
       spc$data$H                  <- NULL

       if(spc$target.name == "H2O"){
           spc$target.name             <- "Water, Liquid"
       }
    }

    return(spc)
}

