sc_time_function=function(
    mod,
    min,
    max,
    length=100,
    cross_section=NULL,
    fill=c(NA,NA),
    alpha=0.3,
    linetype=1,
    linewidth=0.2,
    v_linetype=1,
    v_linewidth=0.2,
    xlab="Follow Up",
    ylab="Outcome",
    main="Plot of Outcome over Time",
    ylim=c(-2,2),
    digits=3,
    center_value=6
){
  
  
  t=seq(min,max,length=length)
  
  
  g1=g2=NULL
  for(e in t){
    g1=cbind(g1,c(1,0,e,0))
    g2=cbind(g2,c(1,1,e,e))
  }
  
  pe1=t(g1)%*%as.matrix(coefficients(summary(mod))[,1])
  pe2=t(g2)%*%as.matrix(coefficients(summary(mod))[,1])
  
  se1=sqrt(diag(t(g1)%*%vcov(mod)%*%g1))
  se2=sqrt(diag(t(g2)%*%vcov(mod)%*%g2))
  
  ll1=pe1+se1*qt(0.025,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  ll2=pe2+se2*qt(0.025,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  
  ul1=pe1+se1*qt(0.975,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  ul2=pe2+se2*qt(0.975,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  
  gg_inf=data.frame(
    m=c(pe1,pe2),
    se=c(se1,se2),
    ll=c(ll1,ll2),
    ul=c(ul1,ul2),
    time=c(t+center_value,t+center_value),
    Group=c(rep("Comparison group",length),rep("TB survivors",length))
  )
  
  
  L2=cbind(c(0,0,1,0),c(0,0,1,1),c(0,0,0,1))
  
  slp=t(L2)%*%as.matrix(coefficients(summary(mod))[,1])
  
  se_slp=sqrt(diag(t(L2)%*%vcov(mod)%*%L2))
  
  ll_slp=slp+se_slp*qt(0.025,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  ul_slp=slp+se_slp*qt(0.975,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  
  p_slp=ifelse(2*pt(-abs(slp)/se_slp,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))<0.001,"< 0.001",round(2*pt(-abs(slp)/se_slp,length(summary(mod)$residuals)-nrow(coefficients(summary(mod)))),3))
  
  blurb=data.frame("Estimates of Changes Over Time")
  rownames(blurb)=""
  colnames(blurb)=""
  
  print(blurb)
  
  out_tab=data.frame(
    PE=round(slp,digits),
    SE=round(se_slp,digits),
    LL=round(ll_slp,digits),
    UL=round(ul_slp,digits),
    p=p_slp
  )
  
  
  
  g=ggplot(data=gg_inf,aes(x=time,y=m,group=Group))+
    geom_line(linetype=linetype,linewidth=linewidth)+
    geom_ribbon(aes(ymin=ll,ymax=ul,fill=Group),alpha=alpha)+
    scale_fill_manual(values=fill)+
    ylim(ylim)+
    xlab(xlab)+
    ylab(ylab)+
    geom_vline(xintercept=cross_section+center_value,linetype=v_linetype,linewidth=v_linewidth)+
    ggtitle(main)
  print(g)
  
  rownames(out_tab)=c("Comparison group","TB survivors","Difference in Slopes")
  print(out_tab)
  
  
  if(is.null(cross_section)==F){
    L=cbind(c(1,0,cross_section,0),c(1,1,cross_section,cross_section))
    L=cbind(L,L[,1]-L[,2])
    delt=t(L)%*%as.matrix(coefficients(summary(mod))[,1])
    se_delt=sqrt(diag(t(L)%*%vcov(mod)%*%L))
    ll_delt=delt+se_delt*qt(0.025,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
    ul_delt=delt+se_delt*qt(0.975,length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
    
    out_tab2=data.frame(
      PE=round(delt,digits),
      SE=round(se_delt,digits),
      LL=round(ll_delt,digits),
      UL=round(ul_delt,digits)
    )
    
    blurb=data.frame(paste("Difference in Averages at ",cross_section+center_value," Month Follow Up",sep=""))
    rownames(blurb)=""
    colnames(blurb)=""
    print(blurb)
    
    rownames(out_tab2)=c("Comparison group","TB survivors","Difference in Averages")
    print(out_tab2)
    
    blurb=data.frame(paste("Test of difference in averages at ",cross_section+center_value," Month follow up: p ",
                           ifelse(2*pt(-abs(delt[3])/se_delt[3],length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))<0.001,"< 0.001",paste("= ",round(2*pt(-abs(delt[3])/se_delt[3],length(summary(mod)$residuals)-nrow(coefficients(summary(mod)))),3),sep="")),sep=""))
    rownames(blurb)=""
    colnames(blurb)=""
    
    print(blurb)
    
  }
  
  blurb=data.frame("Note: PE = Point Estimate; SE = Standard Error; LL = Lower Limit; UL = Upper Limit; p is the p value testing the null hypothesis of zero (i.e., no change over time, or no difference in changes over time)")
  rownames(blurb)=""
  colnames(blurb)=""
  
  print(noquote(blurb))
  
  out=list()
  
  out$slopes=c(slp)
  out$slopes_se=c(se_slp)
  out$slopes_ll=c(ll_slp)
  out$slopes_ul=c(ul_slp)
  out$slopes_t=c(slp/se_slp)
  out$df=length(summary(mod)$residuals)-nrow(coefficients(summary(mod)))
  out$slopes_p=c(2*pt(-abs(slp)/se_slp,length(summary(mod)$residuals)-nrow(coefficients(summary(mod)))))
  
  if(is.null(cross_section)==F){
    out$difference=c(delt)
    out$difference_se=c(se_delt)
    out$difference_ll=c(ll_delt)
    out$difference_ul=c(ul_delt)
    out$difference_t=delt[3]/se_delt[3]
    out$p_value=pt(-abs(delt[3])/se_delt[3],length(summary(mod)$residuals)-nrow(coefficients(summary(mod))))
  }
  
  out$chart_data=gg_inf
  out$chart_code="ggplot(data=gg_inf,aes(x=time,y=m,group=Group))+
  geom_line(linetype=linetype,linewidth=linewidth)+
  geom_ribbon(aes(ymin=ll,ymax=ul,fill=Group),alpha=alpha)+
  scale_fill_manual(values=fill)+
  ylim(ylim)+
  xlab(xlab)+
  ylab(ylab)+
  geom_vline(xintercept=cross_section)+
  ggtitle(main)"
  
  return(out)
}
