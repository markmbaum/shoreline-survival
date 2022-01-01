from os.path import join, isfile
from pandas import read_csv
import matplotlib.pyplot as plt
import matplotx
from seaborn import *
from numpy import *

#------------------------------------------------------------------------------
# INPUT

#csv file with simulation results
fnsim = join('..', '..', 'data', 'sims', 'mapped.csv')

#default overlap (Delta) value
overlap = 50

#directory to save plots in
dirsave = join('..', '..', 'plots', 'results')

#whether or not to save
save = True

#whether to force saving over existing files
overwrite = True

#dots per inch to save with
dpi = 500 #high rez!

#------------------------------------------------------------------------------
#FUNCTIONS

def saveclose(fig, fn):
    path = join(dirsave, fn + '.png')
    if overwrite or not isfile(path):
        fig.savefig(path, dpi=dpi)
        print("figure saved: %s" % path)
        plt.close(fig)

def gyaaxis(ax):
    ax.set_xlabel('Time [Ga]')
    ax.invert_xaxis()

def jitter(v, scale, vmin=None, vmax=None):
    r = scale*(max(v) - min(v))*(2*random.rand(len(v)) - 1)/2
    vj = v + r
    if vmin is not None:
        vj[vj < vmin] = vmin
    if vmax is not None:
        vj[vj > vmax] = vmax
    return vj

def jitterscatter(ax, x, y,
                  xscale=0.05,
                  yscale=0.0,
                  color="C0",
                  size=20,
                  line=True,
                  xmax=None,
                  xmin=None,
                  ymax=None,
                  ymin=None,
                  label=None):
    ax.scatter(
        jitter(x, xscale, xmin, xmax),
        jitter(y, yscale, ymin, ymax),
        color=color,
        s=size,
        edgecolors='k',
        label=label
    )
    if line:
        xl = sort(unique(x))
        yl = [mean(y[x == u]) for u in xl]
        ax.plot(
            xl,
            yl,
            color=color,
            zorder=-1,
            alpha=0.5,
            linewidth=3
        )
    return None

def format(ax):
    leg = ax.legend(fontsize=8, framealpha=1)
    leg.set_title('Ejecta Multiple')
    gyaaxis(ax)
    ax.set_ylabel(None)
    ax.grid(True, axis='y', zorder=-100)
    ax.grid(False, axis='x')
    fig.tight_layout()
    return None

#------------------------------------------------------------------------------
# MAIN

#read the table
df = read_csv(fnsim)

#convert segment columns to kilometers
for col in df.columns:
    if 'seg' in col:
        df[col] /= 1e3

#choose values of overlap and theta for some plots
dg = df[df.overlap == overlap].drop('overlap', axis=1)

#money plot, fraction surviving over time
fig, ax = plt.subplots(1, 1, figsize=(4,3))
G = dg.groupby('re')
colors = plt.cm.cool(linspace(0, 1, len(G.groups)))
for i,k in enumerate(G.groups):
    g = G.get_group(k)
    jitterscatter(
        ax,
        g.t,
        1 - g.f,
        ymax=1,
        ymin=0,
        color=colors[i],
        label=g.re.iat[0]
    )
ax.set_title('Surviving Fraction')
format(ax)
if save:
    saveclose(fig, 'mapped_fraction_surviving')

#maximum segment length over time
fig, ax = plt.subplots(1, 1, figsize=(4,3))
G = dg.groupby('re')
colors = plt.cm.cool(linspace(0, 1, len(G.groups)))
for i,k in enumerate(G.groups):
    g = G.get_group(k)
    jitterscatter(
        ax,
        g.t,
        g.segmax,
        color=colors[i],
        label=g.re.iat[0]
    )
ax.set_title('Maximum Segment Length [km]')
format(ax)
if save:
    saveclose(fig, 'mapped_max_segment')

"""
#check out the number of registered impacts
fig, ax = plt.subplots(1, 1)
stripplot(data=dg, x='t', y='impacts', hue='re', ax=ax)
gyaaxis(ax)
ax.set_ylabel('Number of Impacts')
ax.get_legend().set_title('Ejecta Multiple')
fig.tight_layout()
if save:
    saveclose(fig, 'mapped_impact_count')

#check the effects of overlap threshold
#fg = relplot(data=df, x='overlap', y='f', hue='t', col='re', kind='line', legend=False)
#fg.figure.tight_layout()
#if save:
#    saveclose(fg.figure, 'mapped_overlap_effect')
"""
plt.show()