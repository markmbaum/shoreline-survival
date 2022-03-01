from os.path import join, isfile
from pandas import read_csv
import matplotlib
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
    ax.set_xlabel('Age [Ga]')
    ax.invert_xaxis()
    return None

slice_re = lambda df, re: df[[x in re for x in df.re]]

def jitter(v, scale, vmin=None, vmax=None):
    r = scale*(max(v) - min(v))*linspace(-1, 1, len(v)) #(2*random.rand(len(v)) - 1)/2
    random.shuffle(r)
    vj = v + r
    if vmin is not None:
        vj[vj < vmin] = vmin
    if vmax is not None:
        vj[vj > vmax] = vmax
    return vj

def jitterscatter(ax, x, y,
                  xscale=0.02,
                  color="C0",
                  size=20,
                  line=True,
                  xmax=None,
                  xmin=None,
                  label=None,
                  zorder=1):
    ax.scatter(
        jitter(x, xscale, xmin, xmax),
        y,
        color=color,
        s=size,
        edgecolors='k',
        label=label,
        zorder=zorder
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
            alpha=0.2,
            linewidth=3
        )
    return None

istens = lambda x: round(10*x) == 10*x

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

fig, ax = plt.subplots(1, figsize=(4,3))
lineplot(
    data=dgiso,
    x='t',
    y='survived',
    hue='re',
    ax=ax,
    ci='sd',
    palette=cmap,
    legend=False
)
G = dgmap.groupby('re')
colors = plt.cm.cool(linspace(0, 1, len(G.groups)))
for i,k in enumerate(G.groups):
    g = G.get_group(k)
    jitterscatter(
        ax,
        g.t,
        g.survived,
        color=colors[i],
        line=False,
        zorder=3
    )
gyaaxis(ax)
ax.set_ylim(0, 1)
ax.set_axisbelow(True)
ax.set_title("Shoreline Survival Fraction")
ax.set_ylabel(None)
ax.set_yticks(linspace(0, 1, 6))
ax.set_yticks(linspace(0.1, 0.9, 5), minor=True)
cb = plt.colorbar(
    plt.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(
            dgiso.re.min(),
            dgiso.re.max()
        ),
        cmap=cmap
    ),
    ax=ax
)
cb.set_label(
    "Ejecta Multiple",
    rotation=270,
    va='top'
)
cb.set_ticks([1, 2])
fig.tight_layout()
saveclose(fig, 'survival_fraction')

#--------------------------------------
# characterizing segments and gaps

dh = dgiso[dgiso.t.apply(istens)]
dh = slice_re(dh, (1, 1.5, 2))
cols = ['segmax', 'gapmax']
titles = [
    'Maximum Segment Length [km]',
    'Maximum Gap Length [km]'
]
fns = [
    'segs',
    'gaps'
]
for i in range(2):
    fig, ax = plt.subplots(1, 1, figsize=(4,2.5))
    pointplot(
        data=dh,
        x='t',
        y=cols[i],
        hue='re',
        palette='cool',
        ci='sd',
        errwidth=1.5,
        #dodge=0.2,
        ax=ax
    )
    gyaaxis(ax)
    ax.set_title(titles[i])
    ax.set_ylabel(None)
    leg = ax.get_legend()
    leg.set_title("Ejecta Multiple")
    fig.tight_layout()
    saveclose(fig, fns[i])

#--------------------------------------
# chech the effect of latitude

dh = dfiso[dfiso.overlap == overlap]
dh = dh[dh.t.apply(istens)]
dh = slice_re(dh, (1, 1.5, 2))
g = catplot(
    data=dh,
    x='theta',
    y='survived',
    col='re',
    hue='t',
    ci='sd',
    kind='point'
)
fig = plt.gcf()
axs = fig.axes
xtl = [('$\pi/%d$' % (round(pi/x))) for x in sort(dh.theta.unique())]
for ax,re in zip(axs, (1,1.5,2)):
    ax.set_xlabel(r"$\theta$")
    ax.set_title("Ejecta Multiple = " + str(re))
    ax.set_xticklabels(xtl)
axs[0].set_ylabel("Shoreline Survival Fraction")
g.legend.set_title("Age [Ga]")
saveclose(fig, 'theta_check')


plt.show()