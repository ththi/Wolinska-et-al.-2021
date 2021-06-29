library("car")

###################################
#Correlation growth-microbial load#
###################################
data_FW_load=read.table("FW_microbial_load.txt", h=T)
data_FW_load=data_FW_load[-which(data_FW_load$genotype=="bri301"),]

### change genotype names from old version

data_FW_load[,"genotype"]=gsub("apex1","hub1",data_FW_load[,"genotype"])
data_FW_load[,"genotype"]=gsub("apex2","apex",data_FW_load[,"genotype"])
data_FW_load[,"genotype"]=gsub("apex3","hub2",data_FW_load[,"genotype"])

#####BACTERIA#####
mod=lm(FW_ratio_mean  ~log(X16S_load_mean), data=data_FW_load)
#qqnorm(residuals(mod))
#qqline(residuals(mod))
#shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p1 <- ggplot(data_FW_load, aes(x=log(X16S_load_mean), y= FW_ratio_mean)) +labs(x="Mean B abundance (log)", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.4028, R2=-0.01834 ")
#p1

ggsave( "fig_3_D.pdf", p1)

#####FUNGI#####
mod=lm(FW_ratio_mean  ~log(ITS1_load_mean), data=data_FW_load)
#qqnorm(residuals(mod))
#qqline(residuals(mod))
#shapiro.test(residuals(mod))
Anova(mod)
summary(mod)

p2 <- ggplot(data_FW_load, aes(x=log(ITS1_load_mean), y= FW_ratio_mean)) +labs(x="Mean B abundance (log)", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.005374, R2=0.4196")
#p2
ggsave( "fig_3_E.pdf", p2)
#####OOMYCETES#####
mod=lm(FW_ratio_mean  ~log(oITS1_load_mean), data=data_FW_load)
#qqnorm(residuals(mod))
#qqline(residuals(mod))
#shapiro.test(residuals(mod))
Anova(mod)
summary(mod)

p3 <- ggplot(data_FW_load, aes(x=log(oITS1_load_mean), y= FW_ratio_mean)) +labs(x="Mean O abundance (log)", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.3435, R2=-0.002409")
#p3
ggsave( "fig_3_F.pdf", p3)