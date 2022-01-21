from os.path import join, isfile
from pandas import read_csv
import matplotlib.pyplot as plt
from seaborn import *
from numpy import *
from scipy.interpolate import splrep, splev

#------------------------------------------------------------------------------
# INPUT

#csv files with simulation results
dirsims = join('..', '..', 'data', 'sims')
fniso = join(dirsims, 'isolatitude.csv')
fnmap = join(dirsims, 'mapped.csv')

#default overlap (Delta) value
overlap = 50

#default isolatitude value
theta = 1.0472

#youngest ages to show
agemin = 3.6 #[Ga]

#directory to save plots in
dirsave = join('..', '..', 'plots', 'results')

#whether or not to save
save = True

#whether to force saving over existing files
overwrite = True

#dots per inch to save with
dpi = 400 #high rez!

#colormap to use
cmap = plt.cm.cool

#------------------------------------------------------------------------------
#FUNCTIONS

def saveclose(fig, fn):
    path = join(dirsave, fn + '.png')
    if overwrite or not isfile(path):
        fig.savefig(path, dpi=dpi)
        print("figure saved: %s" % path)
        plt.close(fig)
    return None

def gyaaxis(ax):
    ax.set_xlabel('Time [Ga]')
    ax.invert_xaxis()
    return None

slice_re = lambda df, re: df[[x in re for x in df.re]]

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
        xl = array(sort(unique(x)))
        yl = array([mean(y[x == u]) for u in xl])
        S = splrep(xl, yl)
        xs = linspace(xl.min(), xl.max(), 1000)
        ys = splev(xs, S)
        ax.plot(
            xs,
            ys,
            color=color,
            zorder=-1,
            alpha=0.5,
            linewidth=3
        )
    return None

#------------------------------------------------------------------------------
# MAIN

#read the tables
dfiso = read_csv(fniso)
dfmap = read_csv(fnmap)

#prepare results for plotting
dfiso = dfiso[dfiso.t >= agemin]
dfmap = dfmap[dfmap.t >= agemin]
for col in dfiso.columns:
    if ('seg' in col) or ('gap' in col):
        dfiso.loc[:,col] /= 1e3
        dfmap.loc[:,col] /= 1e3

#filter values for a single slice
dgiso = dfiso[(dfiso.theta == theta) & (dfiso.overlap == overlap)].drop(['theta','overlap'], axis=1)
dgmap = dfmap[dfmap.overlap == overlap].drop('overlap', axis=1)

#--------------------------------------
#fraction survived

fig, (axa, axb) = plt.subplots(1, 2, figsize=(8,4))
lineplot(
    data=dgiso,
    x='t',
    y='survived',
    hue='re',
    ax=axa,
    ci='sd',
    palette=cmap
)
G = dgmap.groupby('re')
colors = plt.cm.cool(linspace(0, 1, len(G.groups)))
for i,k in enumerate(G.groups):
    g = G.get_group(k)
    jitterscatter(
        axb,
        g.t,
        g.survived,
        ymax=1,
        ymin=0,
        color=colors[i],
        label=g.re.iat[0]
    )
leg = axa.get_legend()
leg.set_title('Ejecta Multiple')
axb.legend()
axb.get_legend().set_title('Ejecta Multiple')
axa.set_ylabel('Shoreline Survival Fraction [-]')
axb.set_ylabel(None)
axa.set_title('Isolatitude')
axb.set_title('Mapped')
axb.set_axisbelow(True)
for (c,ax) in zip(('a', 'b'), (axa,axb)): 
    gyaaxis(ax)
    ax.set_ylim(0, 1)
fig.tight_layout()
saveclose(fig, 'survival_fraction')

#--------------------------------------
#

fig, ((axa, axb), (axc, axd)) = plt.subplots(2, 2, figsize=(8,6))

lineplot(
    data=slice_re(dgiso, (1, 1.5, 2)),
    x='t',
    y='segmax',
    hue='re',
    ax=axa,
    ci='sd',
    palette=cmap
)
axa.set_title('Maximum Segment Length [km]')
leg = axa.get_legend()
leg.set_title('Ejecta Multiple')

lineplot(
    data=dgiso,
    x='t',
    y='segmedian',
    hue='re',
    ax=axb,
    ci='sd',
    palette=cmap
)
axb.set_title('Median Segment Length [km]')
leg = axb.get_legend()
leg.set_title('Ejecta Multiple')

lineplot(
    data=slice_re(dgiso, (1, 1.5, 2)),
    x='t',
    y='gapmax',
    hue='re',
    ax=axc,
    ci='sd',
    palette=cmap,
    legend=False
)
axc.set_title('Maximum Gap Length [km]')

lineplot(
    data=dgiso,
    x='t',
    y='gapmedian',
    hue='re',
    ax=axd,
    ci='sd',
    palette=cmap,
    legend=False
)
axd.set_title('Median Gap Length [km]')

for ax in (axa, axb, axc, axd):
    gyaaxis(ax)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
axc.set_xlabel('Time [Ga]')
axd.set_xlabel('Time [Ga]')
fig.tight_layout()
saveclose(fig, 'segments_gaps')

#--------------------------------------

plt.show()