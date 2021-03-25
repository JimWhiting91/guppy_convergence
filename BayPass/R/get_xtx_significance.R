require(corrplot) ; require(ape) ; require(data.table)
#source the baypass R functions (check PATH)

setwd("/gpfs/ts0/home/jw962/HP_LP_trials/baypass")
source("R/baypass_utils.R")

# Read in the SNP info
SNP.info<-data.frame(fread("data/five_aside_STAR_snps"))

args<-commandArgs(TRUE)
river<-as.character(args[1])

# upload estimate of omega
omega=as.matrix(read.table(paste0("data/five_aside_STAR_",river,"_CovMatrix.txt")))
pop.names=c("HP","LP")
dimnames(omega)=list(pop.names,pop.names)

#Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
main=expression("Correlation map based on"~hat(Omega)))

#Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",
main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#Estimates of the XtX differentiation measures
anacore.snp.res=data.frame(rbindlist(lapply(1:16,function(x){
	
	snps<-as.vector(SNP.info[seq(x,nrow(SNP.info),16),])
	
	tmp<-data.frame(fread(paste0("outputs/",river,"_five_aside_STAR_sub",x,"_summary_pi_xtx.out")))
	tmp$snp<-snps
	tmp$order<-seq(x,nrow(SNP.info),16)
	return(tmp)
	})))

anacore.snp.res<-anacore.snp.res[order(anacore.snp.res$order),]
	
#plot(anacore.snp.res$M_XtX)

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
pi.beta.coef_1<-mean(unlist(lapply(1:16,function(x){return(read.table(paste0("outputs/",river,"_five_aside_STAR_sub",x,"_summary_beta_params.out"),h=T)$Mean[1])})))
pi.beta.coef_2<-mean(unlist(lapply(1:16,function(x){return(read.table(paste0("outputs/",river,"_five_aside_STAR_sub",x,"_summary_beta_params.out"),h=T)$Mean[2])})))
pi.beta.coef<-c(pi.beta.coef_1,pi.beta.coef_2)

#upload the original data to obtain total allele count
OG.data<-geno2YN(paste0("data/five_aside_STAR_",river,".geno"))

#Create the POD
setwd("/gpfs/ts0/home/jw962/HP_LP_trials/baypass/outputs")
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=10000,sample.size=OG.data$NN,beta.pi=pi.beta.coef,pi.maf=0,suffix=paste0(river,"_five_aside_STAR_10k_sim"))

#########################################
# Run command line of simulated data, takes ~ 5 mins for 10k SNPs over 16 cores
system(command=paste0("baypass -npop 2 -gfile G.",river,"_five_aside_STAR_10k_sim -outprefix ",river,"_five_aside_STAR_10k_sim -nthreads 16"),intern=F,ignore.stdout=T,ignore.stderr=T,wait=T)
#########################################
# Return
setwd("/gpfs/ts0/home/jw962/HP_LP_trials/baypass")
#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table(paste0("outputs/",river,"_five_aside_STAR_10k_sim_mat_omega.out")))
plot(pod.omega,omega) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega)
#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table(paste0("outputs/",river,"_five_aside_STAR_10k_sim_summary_beta_params.out"),h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)

#######################################################
# XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table(paste0("outputs/",river,"_five_aside_STAR_10k_sim_summary_pi_xtx.out"),h=T)$M_XtX
#compute the 1% threshold
pod.thresh99=quantile(pod.xtx,probs=0.99)
pod.thresh999=quantile(pod.xtx,probs=0.999)

#add the thresh to the actual XtX plot
plot(anacore.snp.res$M_XtX)
abline(h=pod.thresh99,lty=2)

####################################################
# Highlight 'outlier' SNPs
####################################################
outlier_SNPs_99<-anacore.snp.res[anacore.snp.res$M_XtX > pod.thresh99,]
outlier_SNPs_999<-anacore.snp.res[anacore.snp.res$M_XtX > pod.thresh999,]

# Add in SNP info
outlier_SNPs_99<-cbind(outlier_SNPs_99,SNP.info[outlier_SNPs_99$MRK,])
outlier_SNPs_999<-cbind(outlier_SNPs_999,SNP.info[outlier_SNPs_999$MRK,])

# Write to output
write.table(outlier_SNPs_99,paste0("outputs/",river,"_five_aside_STAR_outlier_SNPs_q99.xtx"),row.names=F,sep="\t",quote=F)
write.table(outlier_SNPs_999,paste0("outputs/",river,"_five_aside_STAR_outlier_SNPs_q999.xtx"),row.names=F,sep="\t",quote=F)

####################################################
# Highlight 'outlier' Windows
####################################################
# First separate snp out into chr and loc

# Build windows
window_size<-50000

chrs<-unique(dd$chr)
chrs<-chrs[chrs!="x"]

# Set
dd$BP<-as.integer(dd$BP)

baypass_winds<-data.frame(rbindlist(lapply(chrs,function(x){

tmp<-dd[dd$chr == x,]

winds<-seq(0,max(tmp$BP),window_size)
winds2<-winds+window_size

# Go through windows
window_out<-data.frame(rbindlist(lapply(1:length(winds),function(y){

tmp2<-tmp[tmp$BP < winds2[y] &
            tmp$BP > winds[y],]

if(nrow(tmp2) > 0){
SNP_N<-nrow(tmp2)
outlier_N<-nrow(tmp2[tmp2$BF.dB. >= upperQ,])
} else {
SNP_N<-NA
outlier_N<-NA
}

out<-data.frame(chr=x,
                start=winds[y],
                end=winds2[y],
                SNP_N=SNP_N,
                outlier_N=outlier_N)

return(out)
})))

return(window_out)

})))

# Filter duds
baypass_winds<-na.omit(baypass_winds)

# Calculate cutoff
      SNP_sim<-as.data.frame(seq(1,max(baypass_winds$SNP_N),by=1))
      colnames(SNP_sim)<-"SNP_N"

# We need to know what the probability for expected number of SNPs is which = N of outlier SNPs over total SNPs
      p<-sum(baypass_winds$outlier_N)/sum(baypass_winds$SNP_N)

      # Here we calculate binomial expectation for 0.99 quantile
      for (k in 1:max(baypass_winds$SNP_N)){
        SNP_sim$exp[k]<-qbinom (0.999, k, p)
      }

      #Calculate outliers
      baypass_winds$exp<-qbinom(0.999, baypass_winds$SNP_N, p)

      # Write intermediate output to file
      top_outliers<-baypass_winds[baypass_winds$outlier_N > baypass_winds$exp,]

      # Reorder my strength of effect
      top_outliers$residual_N<-top_outliers$outlier_N-top_outliers$exp
      top_outliers<-top_outliers[order(-top_outliers$residual_N),]

      write.table(top_outliers,
            paste0("outputs/",DATASET,"_auxmodel_baypass_BF_windows_",window_size,"_OUTLIERS.txt"),
            row.names=F,quote=F,sep="\t")

