import matplotlib.pyplot as plt
from .base_plotter import Plotter
import numpy as np

class BarPlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested):
        self.nested = nested
        self.fig, self.axs = plt.subplots()
        if nested:
            self.nested_barchart(contributions)
        else:
            self.barchart(contributions)

        self.style()

    def get_plot(self):
        return self.fig, self.axs
        
    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
        
    def barchart(self, contributions):
        self.axs.bar(contributions.keys(), [contribution.n for contribution in contributions.values()],
            yerr=[contribution.s for contribution in contributions.values()], capsize=5, error_kw={'elinewidth': 0.5})

    def nested_barchart(self, contributions):
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

    def style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.set_ylabel(f"Contribution to {self.index_name}")
        self.axs.set_xlabel("Isotope")
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)