# Descriptive reanalysis of published source tables. Run from R/. No causal parameters are fitted.
DATA <- file.path('..','data')
rd <- function(f) read.csv(file.path(DATA,f),stringsAsFactors=FALSE,encoding='UTF-8')
gas <- rd('sennett_final_incubation_gas.csv'); co2 <- rd('liu_2025_panel_h_co2_efflux.csv')
trapz <- function(x,y) {
  stopifnot(length(x)==length(y),length(x)>1,all(is.finite(x)),all(is.finite(y)),all(diff(x)>0))
  sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
}
window_integral <- function(x,y,a,b) {
  o<-order(x);x<-x[o];y<-y[o];stopifnot(a>=min(x),b<=max(x),a<b)
  xx<-sort(unique(c(a,x[x>a&x<b],b)))
  trapz(xx,approx(x,y,xx,rule=1)$y)
}
M <- do.call(rbind,lapply(split(gas,paste(gas$history,gas$replicate)),function(q) {
  q<-q[order(q$time_h),];stopifnot(all(is.finite(q$N2_cum)),all(diff(q$time_h)>0))
  f<-tail(q$N2_cum,1);i<-which(q$N2_cum>=f/2)[1];stopifnot(f>0,i>1)
  t50<-q$time_h[i-1]+(f/2-q$N2_cum[i-1])*diff(q$time_h[(i-1):i])/diff(q$N2_cum[(i-1):i])
  data.frame(history=q$history[1],replicate=q$replicate[1],N2_half_time_h=t50,
   N2O_time_integral=trapz(q$time_h,q$N2O),N2O_integral_1_24h=window_integral(q$time_h,q$N2O,1,24),
   N2O_peak=max(q$N2O),N2O_final=tail(q$N2O,1),N2_final=f,
   first_time_h=min(q$time_h),last_time_h=max(q$time_h),n_observations=nrow(q),max_sampling_gap_h=max(diff(q$time_h)))
}))
M$history<-factor(M$history,levels=c('Ox','SA','LA'))
rho<-cor(M$N2_half_time_h,M$N2O_time_integral,method='spearman')
W<-do.call(rbind,lapply(c('Paddy soil','Sandy soil'),function(soil) do.call(rbind,lapply(list(c(.5,6),c(.5,48)),function(w) {
 area<-function(tr) {q<-co2[co2$soil==soil & co2$treatment==tr,];window_integral(q$time,q$value,w[1],w[2])}
 u<-area('Original soil');s<-area('Sterilized soil');t<-area('Soil with TBA')
 data.frame(soil=soil,window=sprintf('%g-%g h',w[1],w[2]),untreated_integral=u,
   sterilized_integral=s,tba_integral=t,sterilized_pct=100*s/u,tba_pct=100*t/u)
}))))
agg<-aggregate(cbind(N2O_time_integral,N2O_integral_1_24h,N2_final,N2_half_time_h)~history,M,mean)
write.csv(M,file.path(DATA,'figure3_derived_metrics.csv'),row.names=FALSE)
write.csv(W,file.path(DATA,'figure3_derived_co2_windows.csv'),row.names=FALSE)
write.csv(agg,file.path(DATA,'figure3_history_summary.csv'),row.names=FALSE)
P<-data.frame(history=c('Ox','SA','LA'),published_Table1_integrated_N2O=c(49.9,592,134))
P<-merge(P,agg[,c('history','N2O_time_integral')],by='history',sort=FALSE)
P$note<-'Own trapezoidal integral over each supplied recorded trajectory; published integration bounds/method not fully reconciled, particularly LA. Not an exact reproduction of Table 1.'
write.csv(P,file.path(DATA,'sennett_integration_definition_comparison.csv'),row.names=FALSE)
print(agg);print(W);cat(sprintf('Descriptive pooled Spearman rho = %.3f\n',rho))
