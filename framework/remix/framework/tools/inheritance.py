from pathlib import Path
from types import FunctionType

from remix.framework.api.dataset import DataPackage


def main_inheritance(
    scr_dir: str,
    base_path: str,
    scen_dir: str,
    roundts: int = 0,
    _logger: FunctionType = lambda x: None,
) -> None:
    """Main Inheritance methods. To be used embedded in GAMS code.

    Parameters
    ----------
    scr_dir : str
        Scratch directory where GAMS writes temporary files.
    base_path : str
        Base path of the data directory.
    scen_dir : str
        Path of the scenario directory. The separator should be the one used in unix systems independent of the
        system: '/'
    roundts : int, optional
        If its set to 0, the time series won't be rounded if they are too wide raising an error, by default 0.
    _logger : FunctionType, optional
        A logging function that takes strings as input, by default lambda x: None.
    """
    scenario = Path(scen_dir).as_posix()
    data_package = DataPackage(base_path, logger=_logger)
    data_package.inherit_scenario(scenario, scr_dir, fileformat="csv", roundts=roundts)
