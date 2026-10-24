# Independent base-R verification, retained alongside the author analysis scripts.
rows<-list()
check<-function(label,value,target,tol=1e-6){ok<-is.finite(value)&&abs(value-target)<=tol;rows[[length(rows)+1]]<<-data.frame(check=label,value=value,target=target,tolerance=tol,passed=ok);if(!ok)stop(label)}
C<-function(t,Q,a,k)sum(Q*a*(-expm1(-k*t)))
half<-function(a,k)uniroot(function(t)C(t,c(50,50),a,k)-sum(50*a)/2,c(0,1e5),tol=1e-10)$root
check('P1 half-time',half(c(1,1),c(1,.01)),3.39871749)
check('P1 proportional access',half(c(.5,.5),c(1,.01)),3.39871749)
check('P1 selective access',half(c(1,.1),c(1,.01)),.796745738)
check('P1 rates divided by 3',half(c(1,1),c(1,.01)/3),10.19615247,1e-5)
A<-C(48,120*c(.2,.8),c(.25,.25),c(1,.01));B<-C(48,120*c(.7,.3),c(.7,.7),c(1,.01))
check('SI soil A',A,15.1491986);check('SI soil B',B,68.4066585);check('SI ratio B/A',B/A,4.51552985)
x<-uniroot(function(t)20*(-expm1(-t))-60*(-expm1(-.02*t)),c(1e-5,100),tol=1e-10)$root
check('P2 crossing hours',x,20.2732554)
escape<-function(L){cuts<-sort(unique(c(0,L,12,400)));X0<-0;E<-0;N2<-0
 for(j in seq_len(length(cuts)-1)){a<-cuts[j];b<-cuts[j+1];p<-if(a<12)1 else 0;k<-if(a>=L)1 else 0;r<-.1+k
 f<-function(t)X0*exp(-r*(t-a))+p/r*(-expm1(-r*(t-a)));area<-integrate(f,a,b,rel.tol=1e-11)$value
 E<-E+.1*area;N2<-N2+k*area;X0<-f(b)}
 c(E=E,N2=N2,balance=E+N2+X0)}
e0<-escape(0);e4<-escape(4);e6<-escape(6);e8<-escape(8)
check('P3 escape L0',e0['E'],1.09090909);check('P3 escape L4',e4['E'],1.73018224);check('P3 escape L6',e6['E'],2.44374215)
check('P3 escape L8 rounded',e8['E'],3.36,.005);check('P3 six-hour ratio',e6['E']/e0['E'],2.24009697)
check('P3 N mass balance',e6['balance'],12,1e-9);check('P3 excess L4',e4['E']-e0['E'],.639273146)
check('P3 short-lag approximation',.1*4^2/(2*1.1),.727272727)
check('P3 approximation error percent',100*(.727272727/(e4['E']-e0['E'])-1),13.7655683,1e-5)
check('P3 never-return excess',12-12*.1/1.1,10.9090909)
check('P4 B at R1',1/expm1(1),.581976707);check('P4 post-event at R1',1+1/expm1(1),1.581976707)
check('P4 B threshold',1/expm1(log(2)),1,1e-10)
write.csv(do.call(rbind,rows),file.path('..','data','independent_math_results.csv'),row.names=FALSE)
cat(length(rows),'independent mathematical checks passed\n')
