AT.SPC.read <- function(file.name){
     spc.size            <- numeric(1)
     res                 <- .C( "AT_SPC_get_size_from_filename_R",
                                file.name          = as.character(file.name),
                                spc.size           = as.integer(spc.size),
                                PACKAGE            = "libamtrack")

    depth.step           <- integer(res$spc.size)
    depth.g.cm2          <- numeric(res$spc.size)
    E.MeV.u              <- numeric(res$spc.size)
    DE.MeV.u             <- numeric(res$spc.size)
    particle.no          <- integer(res$spc.size)
    fluence.cm2          <- numeric(res$spc.size)
    n.bins.read          <- integer(1)
    
    res                  <- .C( "AT_SPC_read_data_from_filename_R",
                                file.name          = as.character(file.name),
                                n                  = as.integer(spc.size),
                                depth.step         = as.integer(depth.step),
                                depth.g.cm2        = as.single(depth.g.cm2),
                                E.MeV.u            = as.single(E.MeV.u),
                                DE.MeV.u           = as.single(DE.MeV.u),
                                particle.no        = as.integer(particle.no),
                                fluence.cm2        = as.single(fluence.cm2),
                                n.bins.read        = as.integer(n.bins.read),
                                PACKAGE            = "libamtrack")
    
    return( data.frame(  depth.step             = res$depth.step,
                         depth.g.cm2            = res$depth.g.cm2,
                         E.MeV.u                = res$E.MeV.u,
                         DE.MeV.u               = res$DE.MeV.u,
                         particle.no            = res$particle.no,
                         fluence.cm2            = res$fluence.cm2))
}
