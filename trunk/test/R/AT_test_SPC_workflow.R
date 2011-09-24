require(libamtrack)
rm(list = ls())

path.to.files <- "/home/greilich/workspace/eclipse/libamtrack/data/SPC/*0000.spc"

spc.list <- AT.SPC.get.list(path.to.files = path.to.files)

spc <- AT.SPC.get(spc.list = spc.list, energy.MeV.u = 345)

spc.lower <- AT.SPC.read(spc.list$file.name[3], endian = "little")
spc.upper <- AT.SPC.read(spc.list$file.name[4], endian = "little")

ddd.lower <- AT.SPC.tapply(spc = spc.lower$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.interp <- AT.SPC.tapply(spc = spc$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.upper <- AT.SPC.tapply(spc = spc.upper$spc, INDEX = c("depth.g.cm2"), 
                           FUN = AT.total.D.Gy, additional.arguments = list(c("material.no", 
            "AT.material.no.from.material.name('Water, Liquid')", 
            FALSE), c("stopping.power.source.no", "2", FALSE)), 
        names.results = "D.Gy")

ddd.lower$which <- "lower"
ddd.interp$which <- "interp"
ddd.upper$which <- "upper"
ddd <- rbind.data.frame(ddd.lower, ddd.interp, ddd.upper)

require(lattice)

xyplot(D.Gy ~ depth.g.cm2,
       ddd,
       groups = which,
       auto.key = TRUE,
       type = 'o')

### HERE T BECOMES CLEAR THAT THERE IS A PROBLEM WITH 
### INTERPOLATION in DEPTH --> DOUBLE PEAK

spectrum <- AT.SPC.spectrum.at.depth.g.cm2(spc = spc$spc, depth.g.cm2 = 20.0) 

xyplot(fluence.cm2 ~ E.MeV.u,
       spectrum,
       groups = particle.no,
       type = 's')

