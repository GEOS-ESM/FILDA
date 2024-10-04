import cartopy.crs as ccrs
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
from scipy.stats import kde
from scipy import stats
import numpy as np

try:
    # Try importing Theano and PyMC3
    import theano
    import theano.tensor as tensor
    import pymc3 as pm
except ImportError:
    # Fall back to pytensor and pymc if Theano/PyMC3 are not available
    try:
        import pytensor as theano
        import pytensor.tensor as tensor
        import pymc as pm
    except ImportError:
        raise ImportError("Required libraries are not installed.")

def multiFigure(nRow, nCol, **kwargs):
	"""
	Creates a grid layout of figures with customizable options.
	
	Parameters:
	- nRow: Number of rows in the grid.
	- nCol: Number of columns in the grid.
	- proj: Projection method for the plots (default: ccrs.PlateCarree()).
	- nUint: Number of grid units per figure (default: 25).
	- nColGap, nRowGap: Gaps between columns and rows in grid units (default: 5).
	- nPlot: Total number of plots to display (default: nRow * nCol - 1).
	- figsize: Size of the figure (default: (9, 9)).
	- projPos: Positions to apply the projection (default: []).
	- fontsize: Font size for annotations (default: 22).
	- xlabel, ylabel: Label positions (default: 0.1, 0.975).
	- numOn: Whether to display figure numbers (default: False).
	- cord: Coordinate bounds for the plots (default: [90, -90, -180, 180]).
	
	Returns:
	- A tuple containing the figure, list of axes, and the projection method used.
	"""
	
	# Set default values for optional parameters
	proj = kwargs.get('proj', ccrs.PlateCarree())
	nUint = kwargs.get('nUint', 25)
	nColGap = kwargs.get('nColGap', 5)
	nRowGap = kwargs.get('nRowGap', 5)
	nPlot = kwargs.get('nPlot', nRow * nCol - 1)
	figsize = kwargs.get('figsize', (9, 9))
	projPos = kwargs.get('projPos', [])
	fontsize = kwargs.get('fontsize', 22)
	xlabel = kwargs.get('xlabel', 0.1)
	ylabel = kwargs.get('ylabel', 0.975)
	numOn = kwargs.get('numOn', False)
	cord = kwargs.get('cord', [90, -90, -180, 180])
	
	# Calculate grid size
	nRow_grid = nRow * nUint
	nCol_grid = nCol * nUint
	
	# Create figure and grid specification
	fig = plt.figure(figsize=figsize)
	gs = gridspec.GridSpec(nRow_grid, nCol_grid)
	
	axes = []
	
	for i in range(nRow):
		for j in range(nCol):
			nFigure = i * nCol + j
			
			
			if nFigure <= nPlot:
	
				txt_number = '(' + chr(97 + nFigure) + ')'
				if nFigure in projPos:
					ax = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nColGap, j*nUint:(j + 1)*nUint - nRowGap], projection=proj)
					ax.set_extent(cord, crs=proj)
				else:
					ax = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nColGap, j*nUint:(j + 1)*nUint - nRowGap])
	
				if numOn:
					text = ax.annotate(txt_number, xy=(xlabel, ylabel), xycoords='axes fraction', color='k', ha='right', va='top', fontsize=fontsize)
					text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w'), path_effects.Normal()])
				axes.append(ax)
	
	if len(projPos) > 0:
		return fig, axes, proj
	else:
		return fig, axes

def save_figure(fig, saveName, dpi = 300):
	'''
	Function to save a figure object...
	
	'''

	import matplotlib.pyplot as plt
	
	print(' - Saving', saveName)
	
	fig.savefig(saveName, bbox_inches='tight', dpi=dpi)
	plt.close()

def plot_hdi(ax, samples, map_est = None, hdi_prob = 0.95, weights = 1, nbin = 300, title = None, HDI = True, ETI = False, CI = False, return_stats = False, formatter = '{:.2f}'):

    # Check if weights is a scalar (single number)
    if np.isscalar(weights):
        # If weights is a scalar, create an array of that value with the length of samples
        weights = np.full_like(samples, weights, dtype=np.float)
    else:
        # If weights is not a scalar, ensure it is a numpy array
        weights = np.array(weights)

    stats = {}
    bins = np.linspace(np.nanmin(samples), np.nanmax(samples), nbin)
    kernel = kde.gaussian_kde(samples, weights=weights)

    sample_kde = kernel(bins)

    sample_mean = np.nanmean(samples)

    stats['bins'] = bins
    stats['KDE'] = sample_kde
    stats['mean'] = sample_mean

    if HDI:
        samples_hdi = pm.hdi(samples, hdi_prob=hdi_prob)


    if ETI:
        percentile =  100 - hdi_prob*100       
        samples_eti = [np.percentile(samples, percentile/2.), np.percentile(samples, 100-percentile/2.)]
        print(f"{hdi_prob*100}% ETI: ({samples_eti[0]}, {samples_eti[1]})")

        stats['ETI'] = samples_eti

    if CI:
        mean = np.mean(samples)
        std_dev = np.std(samples, ddof=1)  # Sample standard deviation
        n = len(samples)
        confidence_level = 0.95

        t_critical = stats.t.ppf((1 + confidence_level) / 2, df=n-1)
        margin_error = t_critical * (std_dev / np.sqrt(n))

        lower_bound = mean - margin_error
        upper_bound = mean + margin_error

        print(f"95% CI for the mean: ({lower_bound}, {upper_bound})")    
        stats['ETI'] = [lower_bound, upper_bound]

    ax.plot(bins,sample_kde, color = 'k')

    #     ax.hlines(kernel(np.min(bins)), samples_hdi[0], samples_hdi[1], lw=3, color = 'orange')

    if HDI:
        ax.hlines(0, samples_hdi[0], samples_hdi[1], lw=3, color = 'orange')

        ax.vlines(samples_hdi[0], 0, kernel(samples_hdi[0]), lw=1.5, ls = 'dashed', color = 'gray')
        ax.vlines(samples_hdi[1], 0, kernel(samples_hdi[1]), lw=1.5, ls = 'dashed', color = 'gray')

        ax.plot(samples_hdi[0], kernel(samples_hdi[0]), color = 'gray', marker='s', markersize = 5)
        ax.plot(samples_hdi[1], kernel(samples_hdi[1]), color = 'gray', marker='s', markersize = 5)


    #     ax.hlines(np.nanmin(sample_kde), samples_hdi[0], samples_hdi[1], lw=2, color = 'orange')
        ax.text(samples_hdi[0], kernel(samples_hdi[0])*1.1, formatter.format(samples_hdi[0]), color = 'orange', ha='center')
        ax.text(samples_hdi[1], kernel(samples_hdi[1])*1.1, formatter.format(samples_hdi[1]), color = 'orange', ha='center')

        ax.text((samples_hdi[0] + samples_hdi[1])/2, kernel(np.min(bins))*1.2, '{:.0f}'.format(hdi_prob*100) + '% HDI' , ha='center') #, transform=ax.get_yaxis_transform()

    posXY0      = (1, 1)
    posXY_text0 = (-5, -5)
    equations0 = 'Mean=' + formatter.format(sample_mean)

    if map_est != None:
        equations0 = equations0 + '\nMAP=' + formatter.format(map_est)

    ax.annotate(equations0, xy=posXY0, xytext=posXY_text0, va='top', ha='right', xycoords='axes fraction', textcoords='offset points')    



#     ax.vlines(map_est, np.nanmin(sample_kde), np.nanmax(sample_kde), lw=2, color = 'blue')

    #     ax.hlines(0.001, samples_eti[0], samples_eti[1], lw=2, color = 'blue')
    #     ax.text(samples_hdi[0], 0.0005, '{:.2f}'.format(samples_hdi[0]), ha='left')
    #     ax.text(samples_hdi[1], 0.0005, '{:.2f}'.format(samples_hdi[1]), ha='right')    
    #     ax.set_title('Posterior Density of Ts')

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(True, alpha = 0.2, ls = 'dashed')

    if return_stats:
        return ax, stats
    else:
        return ax


def plot_mcmc_fire(trace, map_estimate=dict(t_s=None, t_f=None, f_s=None, f_f=None)):
    
    print(map_estimate)
    
    fig, axes, _ = multiFigure(2, 2, figsize=[8,6])
    
    print(f" - t_s number of the sample: {len(trace['t_s'])}") 
    axes[0] = plot_hdi(axes[0], trace['t_s'], map_est = map_estimate['t_s'], hdi_prob = 0.95, weights=1, nbin = 500, title = None)
    axes[0].set_xlabel('Ts')
    axes[0].set_yticks([])
    
    print(f" - t_f number of the sample: {len(trace['t_f'])}")
    axes[1] = plot_hdi(axes[1], trace['t_f'], map_est = map_estimate['t_f'], hdi_prob = 0.95, weights=1, nbin = 500, title = None)
    axes[1].set_xlabel('Tf')
    axes[1].set_yticks([])
    
    print(f" - f_s number of the sample: {len(trace['f_s'])}")
    axes[2] = plot_hdi(axes[2], trace['f_s'], map_est = map_estimate['f_s'], hdi_prob = 0.95, weights=1, nbin = 500, title = None, formatter = '{:.3f}')
    axes[2].set_xlabel('Fs')
    axes[2].set_yticks([])
    
    print(f" - f_f number of the sample: {len(trace['f_f'])}")
    axes[3] = plot_hdi(axes[3], trace['f_f'], map_est = map_estimate['f_f'], hdi_prob = 0.95, weights=1, nbin = 500, title = None, formatter = '{:.3f}')
    axes[3].set_xlabel('Ff')
    axes[3].set_yticks([])
    return fig

def plot_mcmc_heatflux(trace):

    fig, axes, _ = multiFigure(1, 3, figsize=[12,3])

    plot_hdi(axes[0], trace['hf_f'], map_est = None, hdi_prob = 0.95, weights = 1, nbin = 300, title = None, HDI = True, ETI = False, CI = False, return_stats = False, formatter = '{:.3f}')

    axes[0].set_xlabel('flaming flux (MW $\cdot$ m$^{-2}$)')
    axes[0].set_yticks([])


    plot_hdi(axes[1], trace['hf_s'], map_est = None, hdi_prob = 0.95, weights = 1, nbin = 300, title = None, HDI = True, ETI = False, CI = False, return_stats = False, formatter = '{:.3f}')
    axes[1].set_xlabel('smoldering flux (MW $\cdot$ m$^{-2}$)')
    axes[1].set_yticks([])


    plot_hdi(axes[2], trace['hf_t'], map_est = None, hdi_prob = 0.95, weights = 1, nbin = 300, title = None, HDI = True, ETI = False, CI = False, return_stats = False, formatter = '{:.3f}')

    axes[2].set_xlabel('total heat flux (MW $\cdot$ m$^{-2}$)')
    axes[2].set_yticks([])
    return fig


def set_month_ticks(ax, is_leap_year=False):
	"""
	Set x-axis to show month labels instead of day of year, adjusting for leap years
	"""
	if is_leap_year:
		# Month ticks for a leap year
		month_ticks = [1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367]
		month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan']
		ax.set_xlim(1, 366)
	else:
		# Month ticks for a non-leap year
		month_ticks = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
		month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan']
		ax.set_xlim(1, 365)
	
	ax.set_xticks(month_ticks)
	ax.set_xticklabels(month_labels)

