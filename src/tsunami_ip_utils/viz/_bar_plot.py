"""Functionality related to matplotlib bar plots of contributions to integral indices."""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from ._base_plotter import _Plotter
import numpy as np
from typing import Dict, Union, Tuple
from uncertainties import ufloat
from tsunami_ip_utils.utils import _isotope_reaction_list_to_nested_dict
import matplotlib as mpl

mpl.rcParams['hatch.linewidth'] = 0.25

class _BarPlotter(_Plotter):
    """Class for creating bar plots of contributions to integral indices on a nuclide-wise and nuclide-reaction-wise basis."""
    _fig: Figure
    """The figure object for the plot."""
    _axs: plt.Axes
    """The axes object for the plot."""
    _index_name: str
    """The name of the integral index being plotted."""
    def __init__(self, integral_index_name: str, plot_redundant: bool=False, **kwargs):
        """Initializes a bar plot of the contributions to the given integral index.
        
        Parameters
        ----------
        integral_index_name
            The name of the integral index whose contributions are to be plotted.
        plot_redundant
            Wether to include redundant/irrelevant reactions in the plot. NOTE: this only applies to nested plots, and
            only affects the plot title; it is expected that the provided data is consistent with the flag.
            
        Notes
        -----
        * Redundant reactions are defined as those which are derived from other reactions, e.g. 'total' and 'capture' reactions
          in SCALE.
        * Irrelevant reactions are defined as those which are not directly cross sections (but rather other nuclear data 
          parameters), e.g. 'chi' in SCALE.
        * A flag for including/excluding redundant/irrelevant reactions was provided since, if the user is expecting the
          the contributions to add up nicely, then the redundant reactions should be excluded, and if only cross sections are
          being considered, then the irrelevant reactions should be excluded.
        """
        self._index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def _create_plot(self, contributions: Union[Dict[str, ufloat], Dict[str, Dict[str, ufloat]]], nested: bool):
        """Creates a bar plot of the given contributions to the integral index.
        
        Parameters
        ----------
        contributions
            * If ``nested`` is ``False``, then this should be a dictionary of the form ``{nuclide: contribution}``, where 
              contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index.
            * If ``nested`` is ``True``, then this should be a dictionary of the form ``{nuclide: {reaction: contribution}}``,
              where contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index 
              through the given reaction.
        nested
            Wether the contributions are on a reaction-wise basis or not."""
        self.nested = nested
        self._fig, self._axs = plt.subplots()
        if nested:
            self._nested_barchart(contributions)
        else:
            self._barchart(contributions)

        self._style()

    def _get_plot(self) -> Tuple[Figure, Axes]:
        return self._fig, self._axs
        
    def _add_to_subplot(self, fig, position) -> Figure:
        return fig.add_subplot(position, sharex=self._axs, sharey=self._axs)
        
    def _barchart(self, contributions: Dict[str, ufloat]):
        """Create a bar chart of the contributions to the integral index on a nuclide-wise basis.
        
        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: contribution}``, where contribution is a ``ufloat`` object representing the
            contribution of the nuclide to the integral index."""
        sorted_contributions = sorted(contributions.items(), key=lambda item: abs(item[1].n), reverse=True)
        labels = [item[0] for item in sorted_contributions]
        values = [item[1].n for item in sorted_contributions]
        errors = [item[1].s for item in sorted_contributions]

        self._axs.bar(labels, values, yerr=errors, capsize=5, error_kw={'elinewidth': 0.5})
        self._axs.axhline(0, color='black', linewidth=0.8)  # Add a line at y=0 for clarity

    def _nested_barchart(self, contributions: Dict[str, Dict[str, ufloat]]):
        """Create a bar chart of the contributions to the integral index on a nuclide-reaction-wise basis."""
        # Collect all unique reactions across all nuclides
        all_reactions = set()
        for nuclide_reactions in contributions.values():
            all_reactions.update(nuclide_reactions.keys())
        all_reactions = sorted(all_reactions)

        # Colors for each reaction type
        cmap = plt.get_cmap('rainbow')
        colors = cmap(np.linspace(0, 1, len(all_reactions)))

        # Calculate total values for each nuclide and sort the data
        total_values = {label: sum(contributions[label][r].n for r in contributions[label]) for label in contributions}
        sorted_data = sorted(total_values.items(), key=lambda x: abs(x[1]), reverse=True)
        sorted_labels = [item[0] for item in sorted_data]
        sorted_total_values = [item[1] for item in sorted_data]

        # Variables to hold the bar positions and labels
        indices = range(len(sorted_labels))
        labels = sorted_labels

        # Bottom offset for each stack
        bottoms_pos = np.zeros(len(sorted_labels))
        bottoms_neg = np.zeros(len(sorted_labels))

        for reaction, color in zip(all_reactions, colors):
            values = [(contributions[nuclide].get(reaction, ufloat(0, 0)).n if reaction in contributions[nuclide] else 0) for nuclide in sorted_labels]
            errs = [(contributions[nuclide].get(reaction, ufloat(0, 0)).s if reaction in contributions[nuclide] else 0) for nuclide in sorted_labels]
            pos_values = [max(0, v) for v in values]
            neg_values = [min(0, v) for v in values]
            self._axs.bar(indices, pos_values, label=reaction, bottom=bottoms_pos, color=color,
                        yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            self._axs.bar(indices, neg_values, bottom=bottoms_neg, color=color,
                        yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            bottoms_pos += pos_values
            bottoms_neg += neg_values

        # Adding 'effective' box with dashed border
        for idx, val in zip(indices, sorted_total_values):
            self._axs.bar(idx, abs(val), bottom=0 if val > 0 else val, color='none', edgecolor='black', hatch='///', linewidth=0.25)

        self._axs.set_xticks(indices)
        self._axs.set_xticklabels(labels)
        self._axs.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')

    def _style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self._index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self._index_name}'
        self._axs.set_ylabel(f"Contribution to {self._index_name}")
        self._axs.set_xlabel("Isotope")
        self._axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self._axs.set_title(title_text)
        ticklabels = self._axs.get_xticklabels()
        self._axs.set_xticks(range(len(ticklabels))) # Have to reset these otherwise the next line gives a warning
        self._axs.set_xticklabels(ticklabels, rotation=45, rotation_mode='anchor', ha='right', va='top') 