AT.FLUKA.read.USRBIN.regs <- function(exp.name, number.of.runs, unit, data.source = 'local', vol.cm3 = NULL, density.g.cm3 = NULL)
{
    for (cur.run in 1:number.of.runs){
        #cur.run <- 1
        file.name            <-    paste(    exp.name,                                        # default: local
                                    AT.SGR.leading.zeros(cur.run, 3),
                                    "_fort.",
                                    unit,
                                    sep = "")
        if(data.source == "condor"){                                                        # condor naming in case cleaning hasn't been run
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
        input                <-    scan(file = file.name, what = "character", strip.white = TRUE, sep = "")
        
        # find no. of particles
        no.of.particles    <-    as.numeric(gsub(",", "", input[grep("followed", input) + 1]))

        # find usrbin outputs
        outputs            <-    grep("Region", input) 

        if (cur.run == 1){
            names                <-    gsub(" ", "", input[outputs + 4])
            bins                <-    as.numeric(input[outputs + 10])
            index                <-    c(1, cumsum(bins)[-length(bins)] + 1)

            # build data.frame
            df            <-    data.frame(    idx                    =        1:sum(bins),
                                                    name                =        character(sum(bins)),
                                                    bins                =        numeric(sum(bins)),
                                                    index                =        numeric(sum(bins)),
                                                    bin                    =        numeric(sum(bins)),
                                                    E.dep.GeV                    =        numeric(sum(bins)),
                                                    E2.dep.GeV2                    =        numeric(sum(bins)))

            class(df$name)        <-    "character"

            for (i in 1:length(names)){
                #i <- 2
                ii                                    <-    df$idx >= index[i] & df$idx < (index[i] + bins[i])
                df$name[ii]                <-    as.character(rep(names[i], sum(ii)))
                df$bins[ii]                <-    rep(bins[i], sum(ii))
                df$index[ii]                <-    rep(index[i], sum(ii))
                df$bin[ii]                <-    1:sum(ii)
                tmp                                    <-    as.numeric(input[outputs[i] + 66 + (1:bins[i]) - 1])
                df$E.dep.GeV                        <-    tmp
                df$E2.dep.GeV2                        <-    tmp^2
            }
        }else{
            for (i in 1:length(names)){
                #i <- 1
                ii                                    <-    df$idx >= index[i] & df$idx < (index[i] + bins[i])
                tmp                                    <-    as.numeric(input[outputs[i] + 66 + (1:bins[i]) - 1])
                df$E.dep.GeV                        <-    df$E.dep.GeV + tmp
                df$E2.dep.GeV2                        <-    df$E2.dep.GeV2 + tmp^2
            }
        }
    }

    if(!(is.null(vol.cm3) | is.null(density.g.cm3))){
        # If no value is given for volume or density, skip dose computation
        # If single value for volume / density is given, apply to all regions the same value
        if(length(vol.cm3) == 1){
            df$vol.cm3           <- rep(vol.cm3, nrow(df))
        }else{
            df$vol.cm3           <- vol.cm3
        }
        if(length(density.g.cm3) == 1){
            df$density.g.cm3     <- rep(density.g.cm3, nrow(df))
        }else{
            df$density.g.cm3     <- density.g.cm3
        }
        
        E.dep.GeV.cm3        <- df$E.dep.GeV / (number.of.runs * df$vol.cm3)
        stdev.E.dep.GeV.cm3  <- sqrt(df$E2.dep.GeV2 / (number.of.runs * df$vol.cm3) - E.dep.GeV.cm3^2)				# valid only for large (>10) numbers of number.of.runs, when 1/N ~ 1/(N-1) for estimating the stdev
        sterr.E.dep.GeV.cm3  <- stdev.E.dep.GeV.cm3 / sqrt(number.of.runs)
        df$D.Gy              <- E.dep.GeV.cm3 * 1.602176462e-7 / df$density.g.cm3
        df$sterr.D.Gy        <- sterr.E.dep.GeV.cm3 * 1.602176462e-7 / df$density.g.cm3
    }

    df$scoring <- rep("Regions", nrow(df))
    return(df)
}
