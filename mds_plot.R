#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
pe <- args[1]
te <- args[2]
pr <- args[3]
de <- args[4]
ra <- args[5]
rb <- args[6]
ga <- args[7]
gb <- args[8]
ma <- args[9]
mb <- args[10]

if (length(args) < 1) {
   print("", quote=F)
   print("Usage: plot.R potener temp press dens rmsda rmsdb gyra gyrb ramaa ramab", quote=F)
   print("", quote=F)
   quit(save="no")
} else {
engy_plot <- function(pe=potential_file,
		     te=temperature_file,
		     pr=pressure_file,
		     de=density_file) {

   in_base <- gsub(".txt","",pe)
   enerPlotOut <- paste0("grom","_energies.png")

   e <- read.table(pe, h=F)
   colnames(e) <- c("Time_ps","Energy_kjmol_1")
   
   t <- read.table(te, h=F)
   colnames(t) <- c("Time_ps","K")
   t_avg <- mean(t$K)
   
   p <- read.table(pr, h=F)
   colnames(p) <- c("Time_ps","bar")
   #p$bar_10ps <- (p$bar)/10
   
   d <- read.table(de, h=F)
   colnames(d) <- c("Time_ps","den")
   #d$den_10ps <- (d$den)/10

   png(enerPlotOut, height=10, width=10, units = "in", points=10, res=200)
   par(mfrow=c(2,2))
      plot(e$Time_ps, e$Energy_kjmol_1, type="l", 
           main=paste0("Potential Energy: ",""," NVT Equilibration"), 
           xlab="Time (ps)", 
           ylab="Potential Energy (kJ/mol)",
           col="red")
      #dev.off()
      
      plot(t$Time_ps, t$K, type="l", 
           main=paste0("Temperature: ",""," NVT Equilibration"), 
           xlab="Time (ps)", 
           ylab="Temperature (K)",
           col = "black")
      abline(h=t_avg, lty=2)
      
      plot(p$Time_ps, p$bar, type="l",
           main=paste0("Pressure: ",""," NPT Equilibration"),
           xlab="Time (ps)",
           ylab="Pressure (bar)",
           col="black")
      #lines(p$Time_ps, p$bar_10ps, col="red")
      
      plot(d$Time_ps, d$den, type="l",
           main=paste0("Density: ",""," NPT Equilibration"),
           xlab="Time (ps)",
           ylab="Density (kg/m^3)",
           col="black")
      #lines(d$Time_ps, d$den_10ps, col="red")
   dev.off()

}

ansis_plot <- function(rmsda=first_rmsd_file,
		       rmsdb=second_rmsd_file, 
		       gyra=first_gyration_file,
		       gyrb=second_gyration_file,
		       ramaa=first_ramachandran,
		       ramab=second_ramachandran) {

   ra <- read.table(rmsda, h=F, col.names=c("time","rmsd"))
   rb <- read.table(rmsdb, h=F, col.names=c("time","rmsd"))

   ga <- read.table(gyra, h=F, col.names=c("time","Rt","Rx","Ry","Rz"))
   ga$t_ns <- ga$time/1000

   gb <- read.table(gyrb, h=F, col.names=c("time","Rt","Rx","Ry","Rz"))
   gb$t_ns <- gb$time/1000

   ma <- read.table(ramaa, h=F, col.names=c("Phi", "Psi", "Residue"))
   mb <- read.table(ramab, h=F, col.names=c("Phi", "Psi", "Residue"))

   png("grom_analysis.png", height=10, width=7, units = "in", points=10, res=200)
   par(mfrow=c(2,1))
      plot(ra$time, ra$rmsd, 
           main="RMSD",
           col="black",
      	   type="l",
           ylim=c(0.15, 0.5),
      	   xlab="Time (ns)",
      	   ylab="RMSD (nm)")
      lines(rb$time, rb$rmsd,col="red")
      legend("bottomright", legend=c("Wild Type","Mutant"), col=c("black","red"),pch="-")

      plot(ga$t_ns, ga$Rt,
           main="Radius of Gyration",
           xlab="Time (ns)",
           ylab="Rg (nm)",
	   ylim=c(1.7, 2.3),
	   col="black",
	   type="l")
      lines(gb$t_ns, gb$Rt,col="red")
      legend("bottomright", legend=c("Wild Type","Mutant"), col=c("black","red"),pch="-")   
   dev.off()
   
   png("ramachandran.png", height=6, width=6, units="in", points=10, res=200)
      plot(ma$Phi, ma$Psi, main="Ramachandran plot", xlab="Phi", ylab="Psi", pch=20)
      points(mb$Phi, mb$Psi, pch=20, col="red")
      abline(v=0, h=0, lty=2)
      legend("bottomright", legend=c("wild type", "mutant"), col=c("black", "red"), pch=20)
   dev.off()
}

engy_plot(pe=pe,te=te,pr=pr,de=de)
ansis_plot(rmsda=ra, rmsdb=rb, gyra=ga, gyrb=gb, ramaa=ma, ramab=mb)
}
