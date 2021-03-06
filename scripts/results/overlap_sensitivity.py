from os.path import join, isfile
from pandas import read_csv
import matplotlib
import matplotlib.pyplot as plt
from seaborn import *
from numpy import *

#------------------------------------------------------------------------------
# INPUT

#csv file with overlap sensitivity results
fn = join('..', '..', 'data', 'sims', 'isolatitude_overlap_sensitivity.csv')

#directory to save plots in
dirsave = join('..', '..', 'plots', 'results')

#dots per inch to save with
dpi = 500 #high rez!

#------------------------------------------------------------------------------
# MAIN

df = read_csv(fn)
df = df[isin(df.overlap.values, [10, 50, 100, 500, 1000, 5000])]
t = sort(df.t.unique())

fg = relplot(
    data=df,
    x='overlap',
    y='survived',
    hue='re',
    col='t',
    height=3,
    kind='line',
    ci='sd',
    palette='cool'
)
fg.axes[0,0].set_ylabel('Survival Fraction')
fg.legend.set_title('Ejecta\nMultiple')
for i,ax in enumerate(fg.axes[0]):
    ax.set_xlabel('Overlap Threshold [m]')
    ax.set_title(str(t[i]) + ' Ga')
    ax.set_xscale('log')
    #ax.set_yscale('log')
fg.figure.savefig(
    join(
        dirsave,
        'overlap_survival_sensitivity'
    ),
    dpi=dpi
)

df['segmax'] = df.segmax/1e3
fg = relplot(
    data=df,
    x='overlap',
    y='segmax',
    hue='re',
    col='t',
    height=3,
    kind='line',
    ci='sd',
    palette='cool'
)
fg.axes[0,0].set_ylabel('Max Segment Length [km]')
fg.legend.set_title('Ejecta\nMultiple')
for i,ax in enumerate(fg.axes[0]):
    ax.set_xlabel('Overlap Threshold [m]')
    ax.set_title(str(t[i]) + ' Ga')
    ax.set_xscale('log')
    ax.set_yscale('log')
fg.figure.savefig(
    join(
        dirsave,
        'overlap_segmax_sensitivity'
    ),
    dpi=dpi
)
