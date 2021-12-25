from os.path import join, isfile
from pandas import read_csv
import matplotlib.pyplot as plt
from seaborn import *
from numpy import *

#------------------------------------------------------------------------------
# INPUT

#csv file with simulation results
fnsim = join('..', '..', 'data', 'sims', 'isolatitude.csv')

#default overlap (Delta) value
overlap = 50

#default isolatitude value
theta = 1.0472

#directory to save plots in
dirsave = join('..', '..', 'plots', 'results')

#whether or not to save
save = True

#whether to force saving over existing files
overwrite = True

#dots per inch to save with
dpi = 400 #high rez!

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

#------------------------------------------------------------------------------
# MAIN

#read the table
df = read_csv(fnsim)

#convert segment columns to kilometers
for col in df.columns:
    if 'seg' in col:
        df[col] /= 1e3

#choose values of overlap and theta for some plots
dg = df[(df.theta == theta) & (df.overlap == overlap)].drop(['theta','overlap'], axis=1)

#money plot, fraction destroyed over time
fig, ax = plt.subplots(1, 1)
lineplot(data=dg, x='t', y='f', hue='re', ax=ax)
gyaaxis(ax)
ax.set_ylabel('Fraction of Hypothetical Shoreline Destroyed')
ax.get_legend().set_title('Ejecta Multiple')
fig.tight_layout()
if save:
    saveclose(fig, 'isolat_fraction_destroyed')

#maximum segment length over time
fig, ax = plt.subplots(1, 1)
lineplot(data=dg, x='t', y='segmax', hue='re', ax=ax)
gyaaxis(ax)
ax.set_ylabel('Maximum Segment Length [km]')
ax.get_legend().set_title('Ejecta Multiple')
fig.tight_layout()
if save:
    saveclose(fig, 'isolat_max_segment')

#check out the number of registered impacts
fig, ax = plt.subplots(1, 1)
lineplot(data=dg, x='t', y='impacts', hue='re', ax=ax)
gyaaxis(ax)
ax.set_ylabel('Maximum Segment Length [km]')
ax.get_legend().set_title('Ejecta Multiple')
fig.tight_layout()
if save:
    saveclose(fig, 'isolat_impact_count')


#check the effects of line latitude and overlap threshold
#fg = relplot(data=df, x='theta', y='f', hue='t', col='re', row='overlap', kind='line', legend=False)
#fg.figure.tight_layout()
#if save:
#    saveclose(fg.figure, 'isolat_colatitude_effect')
#fg = relplot(data=df, x='overlap', y='f', hue='t', col='re', row='theta', kind='line', legend=False)
#fg.figure.tight_layout()
#if save:
#    saveclose(fg.figure, 'isolat_overlap_effect')


plt.show()