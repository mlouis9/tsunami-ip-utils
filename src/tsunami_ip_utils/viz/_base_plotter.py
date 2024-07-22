from abc import ABC, abstractmethod
from typing import Any

class _Plotter(ABC):
    @abstractmethod
    def _create_plot(self, data, nested) -> None:
        pass

    @abstractmethod
    def _add_to_subplot(self, fig: Any, position: tuple) -> Any:
        """Add the plot (corresponding to the class instance) to the given figure subplot at the given position.
        
        Parameters
        ----------
        fig
            The figure to add the plot to.
        position
            The position of the subplot to add the plot to.
            
        Returns
        -------
            The figure with the plot added."""
        pass

    @abstractmethod
    def _get_plot(self) -> Any:
        """Return the 'plot' object.
        
        Returns
        -------
            The plot object."""
        pass

    @abstractmethod
    def _style(self) -> None:
        """Style the plot."""
        pass