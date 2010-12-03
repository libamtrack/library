AT.FLUKA.read.USRBIN.mesh <- function(exp.name, number.of.runs, unit, data.source = 'local', density.g.cm3 = 1.0)
{

    for (cur.run in 1:number.of.runs){
        # DEBUG: cur.run <- 1
        file.name            <-    paste(    exp.name,                
                                    AT.SGR.leading.zeros(cur.run, 3),
                                    "_fort.",
                                    unit,
                                    sep = "")
        if(data.source == "condor"){                                                         # condor naming in case cleaning hasn't been run
            file.name            <-    paste(    exp.name,
                                        "_node",
                                        AT.SGR.leading.zeros(cur.run - 1, 4),
                                        "001_fort.",
                                        unit,
                                        sep = "")
        }
        if(data.source == "condor_cleaned"){
            file.name            <-    paste(    exp.name,                                    # condor
                                        AT.SGR.leading.zeros(cur.run - 1, 5),
                                        "_fort.",
                                        unit,
                                        sep = "")
        }
        input               <-    scan(file = file.name, what = "character", strip.white = T, sep = "")
        
        # find no. of particles
        no.of.particles    <-    as.numeric(gsub(",", "", input[grep("followed", input) + 1]))

        # find usrbin outputs
        outputs            <-    grep("Cartesian", input) 

        # find mesh size
        size               <-    grep("wide", input)
	  sizes.cm           <-    input[size - 2]
        vol.cm3            <-    as.numeric(sizes.cm[1]) * as.numeric(sizes.cm[2]) * as.numeric(sizes.cm[3])

        if (cur.run == 1){
            names               <- gsub(" ", "", input[outputs + 4])
            bins                <- as.numeric(input[outputs + 43])
            index               <- c(1, cumsum(bins)[-length(bins)] + 1)

            # build data.frame
            df    <-    data.frame( idx           = 1:sum(bins),
                                    name          = character(sum(bins)),
                                    bins          = numeric(sum(bins)),
                                    index         = numeric(sum(bins)),
                                    bin           = numeric(sum(bins)),
                                    E.dep.GeV     = numeric(sum(bins)),
                                    E2.dep.GeV2   = numeric(sum(bins)))

            class(df$name)        <-    "character"

            for (i in 1:length(names)){
                #i <- 1
                ii                   <- df$idx >= index[i] & df$idx < (index[i] + bins[i])
                df$name[ii]          <- as.character(rep(names[i], sum(ii)))
                df$bins[ii]          <- rep(bins[i], sum(ii))
                df$index[ii]         <- rep(index[i], sum(ii))
                df$bin[ii]           <- 1:sum(ii)
                tmp                  <- as.numeric(input[outputs[i] + 63 + (1:bins[i]) - 1])
                df$E.dep.GeV         <- tmp
                df$E2.dep.GeV2       <- tmp^2
            }
        }else{
            for (i in 1:length(names)){
                #i <- 1
                ii                   <- df$idx >= index[i] & df$idx < (index[i] + bins[i])
                tmp                  <- as.numeric(input[outputs[i] + 63 + (1:bins[i]) - 1])
                df$E.dep.GeV         <- df$E.dep.GeV + tmp
                df$E2.dep.GeV2       <- df$E2.dep.GeV2 + tmp^2
            }
        }
    }

    df$vol.cm3           <- rep(vol.cm3, nrow(df))
    df$density.g.cm3     <- density.g.cm3
    E.dep.GeV.cm3        <- df$E.dep.GeV / (number.of.runs * df$vol.cm3)
    stdev.E.dep.GeV.cm3  <- sqrt(df$E2.dep.GeV2 / (number.of.runs * df$vol.cm3) - E.dep.GeV.cm3^2)				# valid only for large (>10) numbers of number.of.runs, when 1/N ~ 1/(N-1) for estimating the stdev
    sterr.E.dep.GeV.cm3	 <- stdev.E.dep.GeV.cm3 / sqrt(number.of.runs)
    df$D.Gy              <- E.dep.GeV.cm3 * 1.602176462e-7 / df$density.g.cm3
    df$sterr.D.Gy        <- sterr.E.dep.GeV.cm3 * 1.602176462e-7 / df$density.g.cm3

    df$scoring <- rep("Cartesian mesh", nrow(df))

    return(df)
}
