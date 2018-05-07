f0x, f1x, l0, l1 = fit_result['x']
f0y, _ = analysis(l0, f0x['ini'], 'Sy')
f1y, _ = analysis(l1, f1x['ini'], 'Sy')

##
def plot_hists():
    f, axes = plt.subplots(2,2)
    axes[0,0].hist(f0x['slp'])
    axes[0,0].set_title('no tilt, horiz plane')
    axes[0,1].hist(f1x['slp'])
    axes[0,1].set_title('w/tilt, horiz plane')
    axes[1,0].hist(f0y['slp'])
    axes[1,0].set_title('no tilt, vert plane')
    axes[1,1].hist(f1y['slp'])
    axes[1,1].set_title('w/tilt, vert plane')
    plt.show()
##

def one_one(data, ax=None):
    if ax is None:
        plt.plot(data, data, '--r')
    else:
        ax.plot(data, data, '--r')

def cross_plot():
    def sub(axis, x, y0, y1, title):
        axis.plot(x, y0, '.b', label='no tilt')
        axis.plot(x, y1, '.r', label='w/tilt')
        axis.set_title(title)
        axis.set_xlabel('initial distribution')
        axis.set_ylabel('frequency')
        
    def cross(axis, x, y, title):
        axis.plot(x, y, '.b')
        one_one(x, axis)
        axis.set_title(title)
        axis.set_xlabel('no-tilt frequency')
        axis.set_ylabel('tilted frequency')

    tof = 1e-6
    f, ((ax0, ax1),(ax2, ax3)) = plt.subplots(2,2)
    sub(ax0, f0x['ini'], f0x['slp']/tof, f1x['slp']/tof, 'x-plane')
    cross(ax1, f0x['slp']/tof, f1x['slp']/tof, 'x-plane')
    sub(ax2, f0y['ini'], f0y['slp']/tof, f1y['slp']/tof, 'y-plane')
    cross(ax3, f0y['slp']/tof, f1y['slp']/tof, 'y-plane')
    plt.show()
    
    
