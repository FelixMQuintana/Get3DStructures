from abc import ABC, abstractmethod


class Command(ABC):
    """

    """
    @abstractmethod
    def run(self) -> None:
        """

        :return:
        """
        raise NotImplemented
