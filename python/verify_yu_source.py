"""Independent checks for current Figure 5d and Supplementary Figure 3.
Uses the original bundled Yu et al. workbooks, not just derived CSVs.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import linregress
D=Path(__file__).resolve().parent.parent/"data"
W=D/"primary_source_workbooks"
m4=W/"43247_2026_3307_MOESM4_ESM.xlsx"
m5=W/"43247_2026_3307_MOESM5_ESM.xlsx"
s9=pd.read_excel(m5,sheet_name="Figure S9",header=None)
s5=pd.read_excel(m4,sheet_name="Figure 5",header=None)
s2=pd.read_excel(m4,sheet_name="Figure 2",header=None)
checks=[]
def check(name,value,expected,tol=2e-5):
    a=float(value);b=float(expected);checks.append(dict(check=name,computed=a,expected=b,tolerance=tol,pass_=abs(a-b)<=tol))
d=pd.read_csv(D/"yu_2026_microbial_reduction_eec.csv").set_index("sediment")
for name,row in d.iterrows():
    ss=[i for i,x in s9[1].items() if isinstance(x,str) and x.startswith(name+"-")]
    mos=[i for i,x in s5[0].items() if isinstance(x,str) and x.startswith(name+"-")]
    assert len(ss)==len(mos)==2,(name,ss,mos)
    a,b=s9.iloc[ss[0]],s9.iloc[ss[1]]
    q0,q1=abs(float(a[4])),abs(float(b[4]));r0,r1=float(a[2]),float(b[2])
    val={"EAC_initial":q0,"EAC_end":q1,"EAC_lost":q0-q1,"EAC_lost_sd":math.hypot(a[5],b[5]),
         "EDC_initial":r0,"EDC_end":r1,"EDC_gained":r1-r0,"EDC_gained_sd":math.hypot(a[3],b[3]),
         "clayFeII_pct_initial":s5.iloc[mos[0],1],"clayFeII_pct_end":s5.iloc[mos[1],1]}
    for k,v in val.items():check(f"{name}: {k}",v,row[k])
    check(f"{name}: clay Fe(II) rose >5 percentage points",int(val["clayFeII_pct_end"]-val["clayFeII_pct_initial"]>5),int(row.FeIII_reduced_Mossbauer),0)
def norm(s):return str(s).strip().lower()
f=pd.read_csv(D/"yu_2026_fe_fractions_eac.csv")
for _,r in f.iterrows():
    idx=[i for i,x in s2[0].items() if norm(x)==norm(r.sediment)]
    assert len(idx)==2,(r.sediment,idx)
    a,b=s2.iloc[idx[0]],s2.iloc[idx[1]]
    for col,v in {"Fe_oxyhydroxide_fractions":float(a[1])+float(a[3])+float(a[5]),
                  "Fe_AAO_phyllosilicate_fractions":float(a[7])+float(a[9]),
                  "HCl_FeIII":float(b[7]),"HCl_FeIII_sd":float(b[8]),
                  "EAC":abs(float(b[9])),"EAC_sd":float(b[10])}.items():
        check(f"{r.sediment}: {col}",v,r[col])
for col,R2,p in [("Fe_oxyhydroxide_fractions",.497,.007),("Fe_AAO_phyllosilicate_fractions",.220,.106),("HCl_FeIII",.557,.003)]:
    lr=linregress(f[col],f.EAC);check(f"{col}: R2",lr.rvalue**2,R2,.0006);check(f"{col}: P",lr.pvalue,p,.0006)
check("Mean EAC / HCl Fe(III)",(f.EAC/f.HCl_FeIII).mean(),.93,.005)
out=pd.DataFrame(checks)
out.to_csv(D/"verification_yu_source.csv",index=False)
bad=out[~out.pass_]
print(f"Yu et al.: {len(out)-len(bad)}/{len(out)} checks passed")
if len(bad):
    print(bad.to_string(index=False));raise SystemExit(1)
