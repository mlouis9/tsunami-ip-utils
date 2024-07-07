from abc import ABC, abstractmethod

class Plotter(ABC):
    @abstractmethod
    def create_plot(self, data, nested):
        pass

    @abstractmethod
    def add_to_subplot(self, fig, position):
        pass

    @abstractmethod
    def get_plot(self):
        pass

    @abstractmethod
    def style(self):
        pass