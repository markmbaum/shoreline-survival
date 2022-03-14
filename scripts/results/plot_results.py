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

#default overlap (Delta) values
deltas = [50, 500]

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
dgiso = dfiso[dfiso.theta == theta].drop('theta', axis=1)

#--------------------------------------
#fraction survived

fig, axs = plt.subplots(1, 2, figsize=(7.5,3))
for i,ax,delta in zip(range(2), axs, deltas):
    sliso = dgiso[dgiso.overlap == deltas[i]]
    slmap = dfmap[dfmap.overlap == deltas[i]]
    lineplot(
        data=sliso,
        x='t',
        y='survived',
        hue='re',
        ax=ax,
        ci='sd',
        palette=cmap,
        legend=False
    )
    G = slmap.groupby('re')
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
    ax.set_ylabel(None)
    ax.set_yticks(linspace(0, 1, 6))
    ax.set_yticks(linspace(0.1, 0.9, 5), minor=True)
    ax.set_title(r'$\Delta = %g$ m' % delta)

axs[0].set_ylabel("Shoreline Survival Fraction")
axs[1].set_yticklabels([])
fig.tight_layout()
cb = fig.colorbar(
    plt.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(
            dgiso.re.min(),
            dgiso.re.max()
        ),
        cmap=cmap
    ),
    ax=axs,
    pad=0.02
)
cb.set_label(
    "Ejecta Multiple",
    rotation=270,
    va='top'
)
cb.set_ticks([1, 2])
saveclose(fig, 'survival_fraction')

#--------------------------------------
# characterizing segments and gaps

dh = dgiso[dgiso.t.apply(istens)]
dh = slice_re(dh, (1, 1.5, 2))
cols = ['segmax', 'gapmax']
ylabels = [
    'Maximum Segment Length [km]',
    'Maximum Gap Length [km]'
]
fns = [
    'segs',
    'gaps'
]
for i in range(2):
    fig, axs = plt.subplots(2, 1, figsize=(4,4.5))
    for ax,delta in zip(axs,deltas):
        pointplot(
            data=dh[dh.overlap == delta],
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
        leg = ax.get_legend()
        leg.set_title("Ejecta Multiple")
        ax.set_title(None)
        ax.set_ylabel(None)
        ax.set_title(r'$\Delta = %g$ m' % delta)
    axs[0].set_xlabel(None)
    fig.supylabel(ylabels[i])
    fig.tight_layout()
    saveclose(fig, fns[i])

#--------------------------------------
# check the effect of latitude

dh = dfiso[dfiso.overlap == deltas[0]]
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