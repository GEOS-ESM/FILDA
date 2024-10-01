# import theano.tensor as tensor

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

import numpy as np
from . import MBFPE_MODEL_TENSOR as MT # module for 
from . import MBFPE_VIS as VIS # module for 


def do_eval_uniphasic(est_fire, ts_fire, fss):

	# Define Theano variables for inputs
	t_mean = tensor.scalar('t_mean')
	t_mean.tag.test_value = np.array(0.0, dtype=np.float32)
	
	p_mean = tensor.scalar('p_mean')
	p_mean.tag.test_value = np.array(0.0, dtype=np.float32)
	
	
	uniphasic_radiance = MT.get_band_uniphasic_radiance_tensor(t_mean, p_mean, ts_fire, fss)
	output = uniphasic_radiance.eval({t_mean: est_fire['t_mean'],
									  p_mean: est_fire['p_mean']})

	return output
	


def do_eval_biphasic(est_fire, ts_fire, fss):

	# Define Theano variables for inputs
	t_s = tensor.scalar('t_s')
	t_s.tag.test_value = np.array(0.0, dtype=np.float32)

	t_f = tensor.scalar('t_s')
	t_f.tag.test_value = np.array(0.0, dtype=np.float32)
	
	p_s = tensor.scalar('p_s')
	p_s.tag.test_value = np.array(0.0, dtype=np.float32)
	
	p_f = tensor.scalar('p_s')
	p_f.tag.test_value = np.array(0.0, dtype=np.float32)
	
	
	
	biphasic_radiance = MT.get_band_biphasic_radiance_tensor(t_s, t_f, p_s, p_f, ts_fire, fss)
	output = biphasic_radiance.eval({t_s: est_fire['t_s'],
									 t_f: est_fire['t_f'],
									 p_s: est_fire['p_s'],
									 p_f: est_fire['p_f']})

	return output

	
	
def do_eval_background(est_bg, ts_fire, fss):

	# Define Theano variables for inputs
	t_b = tensor.scalar('t_b')
	t_b.tag.test_value = np.array(0.0, dtype=np.float32)
	
	C = tensor.scalar('C')
	C.tag.test_value = np.array(0.0, dtype=np.float32)
	
	
	uniphasic_radiance = MT.get_band_radiance_BG_tensor(t_b, C, ts_fire, fss)
	output = uniphasic_radiance.eval({t_b: est_bg['t_b'],
									  C  : est_bg['C']})

	return output



def do_eval(nl, ts_bg, ts_bi_fire, ts_uni_fire, bgs, fss, 
		    this_point, data_bg, data_fire, 
		    est_bg, est_fire,
		    savename='eval.png'):

	wl_fire = np.nanmean(fss.lambdas, axis = 1)
	wl_bg   = np.nanmean(bgs.lambdas, axis = 1)
	
	fire_obs = data_fire.raw[this_point]
	bg_obs   = data_bg[this_point][0]
	
	if est_fire['flag_mode'] == 1:
		eval_fire  = do_eval_uniphasic(est_fire, ts_uni_fire, fss)
		eval_background = do_eval_background(est_bg, ts_uni_fire, fss)
		title = f'Uniphasic evaluation'
	else:
		eval_fire  = do_eval_biphasic(est_fire, ts_bi_fire, fss)
		eval_background = do_eval_background(est_bg, ts_bi_fire, fss)
		title = f'Biphasic evaluation'
		
		

	model_sig = eval_background + eval_fire
	
	eval_background_2 = do_eval_background(est_bg, ts_bg, bgs)
	
	
	fig, axes = VIS.multiFigure(1, 1, figsize = (4.5,3))
	
	axes[0].plot(wl_bg, eval_background_2, color = 'C4', 
	             label = 'Modeled background rad')
	axes[0].plot(wl_bg, bg_obs, color = 'C0', 
	             label = 'Obs. background rad')
	
	axes[0].plot(wl_fire, model_sig, color = 'C1', lw=2,
	             label = 'Modeled signal', zorder = 10)
	axes[0].plot(wl_fire, fire_obs, color = 'C3', lw=2, 
	             label = 'Obs. signal', zorder = 10)

	axes[0].set_ylabel('Radiance (W$\cdot$ m$^{-2} \cdot\mu$m$^{-2}$)')
	axes[0].set_xlabel('Wavelength ($\mu m$)')
	axes[0].legend(frameon = False, bbox_to_anchor=(1, 0.7))

	axes[0].set_title(title)

		
	VIS.save_figure(fig, savename)

	
	return
	
	
	
	
	
	
	
	
	
	
	
	