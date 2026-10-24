# Finite-window weak-escape sensitivity: E_W(l)=l*A_W(0)+O(l^2), prescribed source/sink.
source('theory_4R.R')
rows<-list();n<-0L
piece_integral<-function(fun,W,L,T=12) {
 cuts<-sort(unique(c(0,min(L,W),min(T,W),W)));sum(vapply(seq_len(length(cuts)-1),function(i)
 integrate(fun,cuts[i],cuts[i+1],rel.tol=1e-10,abs.tol=1e-11)$value,numeric(1)))
}
for(L in c(0,4,14,30))for(W in c(6,24)) {
 A0<-piece_integral(function(t)x_traj(t,L,l=0),W,L)
 ells<-c(1e-4,5e-5,2.5e-5)
 E<-vapply(ells,function(ell)ell*piece_integral(function(t)x_traj(t,L,l=ell),W,L),numeric(1))
 relerr<-abs(E/(ells*A0)-1)
 passed<-all(is.finite(E))&&all(E>0)&&all(E<=ells*A0+1e-9)&&all(diff(relerr)<0)&&relerr[3]<.001
 n<-n+1L;rows[[n]]<-data.frame(lag_h=L,window_end_h=W,A0=A0,min_leakage=ells[3],
 relative_error_at_min=relerr[3],error_contraction=relerr[3]/relerr[2],ok=passed)
}
D<-do.call(rbind,rows)
# Independent zero-leak construction before/after formation and sink activation.
expected<-c(0,4,12,12*exp(-1))
stopifnot(max(abs(x_traj(c(0,4,12,13),L=12,l=0)-expected))<1e-10)
stopifnot(all(D$ok))
write.csv(D,file.path('..','data','weak_escape_verification_R.csv'),row.names=FALSE)
print(D);cat('8 finite-window cases and zero-leak solution passed. Not a fit to incubation data.\n')
