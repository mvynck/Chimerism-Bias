par(mfrow=c(2,3))
b.max=0.06
b.vec <- seq(b.max, 0.0, -b.max/3)
b.vec
c=seq(0,1,0.0001)
onePlot <- data.frame(bias = c(), obsbias = c(), constellation = c(), HC = c())


# Aa donor, AA recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = (1+c)*(1-2*b)/((1+c )*(1-2*b) + (1-c)*(1+2*b))
  est = 1 - (1 - ref)*2
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="Aa donor\nAA recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-II: Aa/AA", HC = 100*c))
}

# aa donor, AA recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = c*(1-2*b)/(c*(1-2*b) + (1-c)*(1+2*b))
  est = ref
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="aa donor\nAA recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-I: aa/AA", HC = 100*c))
}

# AA donor aa recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = (1-c)*(1-2*b)/((1-c)*(1-2*b)+(c)*(1+2*b))
  est = 1 - ref
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="AA donor\naa recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-I: AA/aa", HC = 100*c))
  
}

# Aa donor aa recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = (1-c)*(1-2*b)/((1-c)*(1-2*b)+(1+c)*(1+2*b))
  est = 1 - ref*2
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="Aa donor\naa recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-II: Aa/aa", HC = 100*c))
}

# aa donor Aa recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = (c)*(1-2*b)/((c)*(1-2*b)+(2-c)*(1+2*b))
  est = ref * 2
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="aa donor\nAa recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-I: aa/Aa", HC = 100*c))
}

# AA donor Aa recipient
for(i in 1:length(b.vec)){
  b = b.vec[i]
  ref = (2-c)*(1-2*b)/((2-c)*(1-2*b)+(c)*(1+2*b))
  est = (1 - ref) * 2
  if(i ==1){
    plot(100*c, 100*(c-est), type = "l", xlab = "True host chimerism (%)", ylab = expression("Observed bias ("*HC[b] - HC[t]*", %)"), main="AA donor\nAa recipient")
  } else {
    lines(100*c, 100*(c-est), lty=i)
  }
  onePlot <- rbind(onePlot, data.frame(bias=b.vec[i], obsbias=100*(c-est), constellation = "type-I: AA/Aa", HC = 100*c))
  
}
onePlot$bias <- as.factor(onePlot$bias)

# combine in one plot
library(ggplot2)
a=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[1],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
b=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[2],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
c=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[3],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
d=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[4],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
e=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[5],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
f=ggplot(onePlot[onePlot$constellation == unique(onePlot$constellation)[6],], aes(x=HC, y=obsbias, group=bias))+
  geom_line(aes(colour=bias))+
  theme_bw()+
  ylab(expression("Observed bias ("*HC[b] - HC[t]*", %)"))+xlab("True host chimerism (%)")+
  scale_colour_manual(values=c("#000000", "#555555", "#999999", "#CCCCCC")) +
  theme(panel.border = element_blank(), legend.position = "none")
pdf("~/Desktop/Figure2.pdf", width = 12, height = 6)
plot_grid(plotlist=list(a,b,c,d,e,f), labels = "AUTO")
dev.off()