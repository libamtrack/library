AT.SPC.read <- function(file.name, mean = c("arithmetic", "geometric")[1], compress = TRUE){
     spc.size            <- numeric(1)
     res                 <- .C( "AT_SPC_get_size_from_filename_R",
                                file.name          = as.character(file.name),
                                spc.size           = as.integer(spc.size),
                                PACKAGE            = "libamtrack")

    n                    <- res$spc.size
    depth.step           <- integer(n)
    depth.g.cm2          <- numeric(n)
    E.MeV.u              <- numeric(n)
    DE.MeV.u             <- numeric(n)
    particle.no          <- integer(n)
    fluence.cm2          <- numeric(n)
    n.bins.read          <- integer(1)
    
    res                  <- .C( "AT_SPC_read_data_from_filename_R",
                                file.name          = as.character(file.name),
                                n                  = as.integer(n),
                                depth.step         = as.integer(depth.step),
                                depth.g.cm2        = as.single(depth.g.cm2),
                                E.MeV.u            = as.single(E.MeV.u),
                                DE.MeV.u           = as.single(DE.MeV.u),
                                particle.no        = as.integer(particle.no),
                                fluence.cm2        = as.single(fluence.cm2),
                                n.bins.read        = as.integer(n.bins.read),
                                PACKAGE            = "libamtrack")
    
    df   <- data.frame(  depth.step             = res$depth.step,
                         depth.g.cm2            = res$depth.g.cm2,
                         E.MeV.u                = res$E.MeV.u,
                         DE.MeV.u               = res$DE.MeV.u,
                         particle.no            = res$particle.no,
                         fluence.cm2            = res$fluence.cm2)

    if (compress == TRUE){
        df  <- df[df$fluence.cm2 != 0,]
    }
}
