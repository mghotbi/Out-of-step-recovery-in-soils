"""Model definitions follow the manuscript: M0 = bulk conditions, pre-challenge (present) function and
challenge forcing; M1 = M0 + internal state; MH = M0 + history; M2 = M0 + both.

Reproduce Figure 6b's illustrative numerical comparison and retain the
exact author-supplied full Figure 6 artwork from the revised manuscript.

This is a proposed experiment, not a fitted model comparison. The lower and
upper bounds in the CSV match the labels in the final BioRender figure.
The unlabelled middle comparison's point is set to -0.17 for the independent
reproduction of panel b; it is schematic and not an empirical estimate.
"""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from style import plt, save
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"figures"
OUT.mkdir(exist_ok=True)
df=pd.read_csv(ROOT/"data"/"fig6_illustrative_intervals.csv")
# Numerical classification of the illustrated intervals.
delta,eps=.10,.05
assert (df.iloc[:2].upper < -delta).all(), "Illustrative state gains must exceed the stated margin"
assert df.iloc[2].lower >= -eps and df.iloc[2].upper <= eps
assert ((df.lower <= df.estimate) & (df.estimate <= df.upper)).all()
fig,ax=plt.subplots(figsize=(5.0,2.25))
fig.subplots_adjust(left=.31,right=.97,bottom=.24,top=.84)
ax.axvspan(-eps,eps,color="#DBDBDC",zorder=0)
ax.axvline(-delta,color="#444D58",ls="--",lw=.8)
ax.axvline(0,color="#888888",lw=.6)
cols=["#CC895F","#A37965","#87529B"]; labels=["−0.22",""," +0.01"]
def fmt(x):
    if x==0:return "0"
    return f"{x:+.2f}".replace("-","−")
for i,row in df.iterrows():
    y=2-i
    ax.errorbar(row.estimate,y,xerr=[[row.estimate-row.lower],[row.upper-row.estimate]],
                fmt="o",ms=5.2,color=cols[i],elinewidth=1.6,capsize=3,zorder=3)
    if labels[i].strip():ax.text(row.estimate,y+.18,labels[i],color=cols[i],ha="center",va="bottom",fontweight="bold")
    ax.text(row.lower,y-.18,fmt(row.lower),ha="center",va="top")
    ax.text(row.upper,y-.18,fmt(row.upper),ha="center",va="top")
ax.set(yticks=[2,1,0],yticklabels=["M1 − M0\n(+ internal state)","M2 − MH\n(+ state, history known)","M2 − M1\n(+ history)"],
       xlim=(-.40,.20),ylim=(-.42,2.58),xlabel="Change in held-out prediction error\n(MAE, log$_{10}$ cumulative N$_2$O-N)")
ax.tick_params(axis="y",length=0)
ax.spines["left"].set_visible(False)
ax.set_xticks(np.arange(-.4,.201,.1))
ax.text(-delta,2.62,"−δ (−0.10)",ha="right",va="bottom")
ax.text(0,2.62,"±ε (±0.05)",ha="center",va="bottom")
save(fig,str(OUT),"Figure_6b_intervals")
# Complete figure includes author-created BioRender panels; preserve the
# exact approved artwork rather than pretending it is fully code-generated.
source=ROOT/"artwork"/"Figure_6_final_author_artwork.png"
shutil.copyfile(source,OUT/"Figure_6_test.png")
print("Figure 6: preserved exact full artwork; Figure 6b interval checks passed.")
