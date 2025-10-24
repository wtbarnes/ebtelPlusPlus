"""
Produce figures for JOSS paper
"""
import astropy.units as u
import matplotlib.pyplot as plt
import pathlib

import ebtelplusplus

total_time = 1.6e4 * u.s
loop_length = 40 * u.Mm

# Setup heating models (all have same event, different partition)
heating_event = ebtelplusplus.models.TriangularHeatingEvent(250*u.s, 200*u.s, 0.05*u.Unit('erg cm-3 s-1'))
equal_heating = ebtelplusplus.models.HeatingModel(events=[heating_event])
electron_heating = ebtelplusplus.models.HeatingModel(partition=1.0,events=[heating_event])
ion_heating = ebtelplusplus.models.HeatingModel(partition=0,events=[heating_event])
single_fluid = ebtelplusplus.models.PhysicsModel(force_single_fluid=True)

# Setup different physics models
expansion = ebtelplusplus.models.PhysicsModel(force_single_fluid=True,
                                              area_ratio_tr_corona=1/3,
                                              area_ratio_0_corona=1/3,
                                              loop_length_ratio_tr_total=0.15)
variable_abundances = ebtelplusplus.models.PhysicsModel(force_single_fluid=True, radiative_loss='variable')

# Setup model combinations
models = {
    'single-fluid': {'heating': equal_heating, 'physics': single_fluid},
    'electron heating': {'heating': electron_heating},
    'ion heating': {'heating': ion_heating},
    'area expansion': {'heating': equal_heating, 'physics': expansion},
    'variable abundance': {'heating': equal_heating, 'physics': variable_abundances},
}

# Build figure
fig,axes = plt.subplot_mosaic(
    [['heating','temperature'],['density', 'phase_space']],
    figsize=(8,4),
    layout='constrained'
)
for model_name, params in models.items():
    result = ebtelplusplus.run(total_time, loop_length, **params)
    l, = axes['temperature'].plot(result.time, result.electron_temperature.to('MK'), label=model_name)
    axes['temperature'].plot(result.time, result.ion_temperature.to('MK'), color=l.get_color(), ls='--')
    axes['density'].plot(result.time, result.density/1e9, color=l.get_color())
    axes['phase_space'].plot(result.electron_temperature.to('K'), result.density, color=l.get_color())
    axes['phase_space'].plot(result.ion_temperature.to('K'), result.density, color=l.get_color(), ls='--')
axes['heating'].plot(result.time, result.heat, color='k')
axes['temperature'].legend(loc=1)
axes['phase_space'].set_xscale('log')
axes['phase_space'].set_yscale('log')
for ax_name in ['heating','temperature','density']:
    axes[ax_name].set_xlim(1e2,total_time.to_value('s'))
    axes[ax_name].set_xlabel('time [s]')
    axes[ax_name].set_xscale('log')
axes['heating'].set_ylabel(r'Heating [erg cm$^{-3}$ s$^{-1}$]')
axes['heating'].set_ylim(1e-7,0.052)
axes['temperature'].set_ylabel(r'Temperature [MK]')
axes['temperature'].set_ylim(0.01,23)
axes['density'].set_ylabel(r'Number Density [$10^9$ cm$^{-3}$]')
axes['density'].set_ylim(1e-3,3)
axes['phase_space'].set_xlabel('Temperature [K]')
axes['phase_space'].set_ylabel(r'Number Density [cm$^{-3}$]')
fig.savefig(pathlib.Path(__file__).parent / 'figure.pdf')
