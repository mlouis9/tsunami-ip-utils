"""Functionality related to matplotlib bar plots of contributions to integral indices."""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from ._base_plotter import _Plotter
import numpy as np
from typing import Dict, Union, Tuple
from uncertainties import ufloat

class _BarPlotter(_Plotter):
    """Class for creating bar plots of contributions to integral indices on a nuclide-wise and nuclide-reaction-wise basis."""
    fig: Figure
    """The figure object for the plot."""
    axs: plt.Axes
    """The axes object for the plot."""
    index_name: str
    """The name of the integral index being plotted."""
    def __init__(self, integral_index_name: str, plot_redundant: bool=False):
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
        self.index_name = integral_index_name
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
        self.fig, self.axs = plt.subplots()
        if nested:
            self._nested_barchart(contributions)
        else:
            self._barchart(contributions)

        self._style()

    def _get_plot(self) -> Tuple[Figure, plt.Axes]:
        return self.fig, self.axs
        
    def _add_to_subplot(self, fig, position) -> Figure:
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
        
    def _barchart(self, contributions: Dict[str, ufloat]):
        """Create a bar chart of the contributions to the integral index on a nuclide-wise basis.
        
        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: contribution}``, where contribution is a ``ufloat`` object representing the
            contribution of the nuclide to the integral index."""
        self.axs.bar(contributions.keys(), [contribution.n for contribution in contributions.values()],
            yerr=[contribution.s for contribution in contributions.values()], capsize=5, error_kw={'elinewidth': 0.5})

    def _nested_barchart(self, contributions):
        """Create a bar chart of the contributions to the integral index on a nuclide-reaction-wise basis.

        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: {reaction: contribution}}``, where contribution is a ``ufloat`` object
            representing the contribution of the nuclide to the integral index through the given reaction."""
        
        # Colors for each reaction type
        num_reactions = len(next(iter(contributions.values())))
        cmap = plt.get_cmap('Set1')
        colors = cmap(np.linspace(0, 1, num_reactions))

        # Variables to hold the bar positions and labels
        indices = range(len(contributions))
        labels = list(contributions.keys())

        # Bottom offset for each stack
        bottoms_pos = [0] * len(contributions)
        bottoms_neg = [0] * len(contributions)

        color_index = 0
        for reaction in next(iter(contributions.values())).keys():
            values = [contributions[nuclide][reaction].n for nuclide in contributions]
            errs = [contributions[nuclide][reaction].s for nuclide in contributions]
            # Stacking positive values
            pos_values = [max(0, v) for v in values]
            neg_values = [min(0, v) for v in values]
            self.axs.bar(indices, pos_values, label=reaction, bottom=bottoms_pos, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            self.axs.bar(indices, neg_values, bottom=bottoms_neg, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            # Update the bottom positions
            bottoms_pos = [bottoms_pos[i] + pos_values[i] for i in range(len(bottoms_pos))]
            bottoms_neg = [bottoms_neg[i] + neg_values[i] for i in range(len(bottoms_neg))]
            color_index += 1

        # Adding 'effective' box with dashed border
        total_values = [sum(contributions[label][r].n for r in contributions[label]) for label in labels]
        for idx, val in zip(indices, total_values):
            self.axs.bar(idx, abs(val), bottom=0 if val > 0 else val, color='none', edgecolor='black', hatch='///', linewidth=0.5)

        self.axs.set_xticks(indices)
        self.axs.set_xticklabels(labels)
        self.axs.legend()

    def _style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.set_ylabel(f"Contribution to {self.index_name}")
        self.axs.set_xlabel("Isotope")
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)