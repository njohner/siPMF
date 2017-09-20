"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This file contains the :class:`CollectiveVariable` class.
"""


class CollectiveVariable():
  """
  This class represents a collective variable (CV or reaction coordinate).
  The name of the CV will give rise to the corresponding fields that are replaced
  in the MD input files. The minimal and maximal values determine the range in which
  the windows can be generated. The number of bins are used for the calculation of the PMF.
  The periodicity of the collective variables is used both for generating new windows and for the
  calculation of the PMF. The step size defines the distance between two neighboring windows.
  Finally the units are just used in the labels of the plots of the PMF.
  """

  def __repr__(self):
    return "CollectiveVariable({0},{1},{2},{3},{4},{5})".format(self.name, self.min_value, self.max_value, self.step_size, self.num_bins, self.periodicity, self.units)

  def __init__(self, cv_name, min_value, max_value, step_size, num_bins, min_spring_constant, max_spring_constant=None, max_shift=None, periodicity=None, units=""):
    """
    :param cv_name: The name of the CV
    :param min_value: The minimal value of the CV
    :param max_value: The maximal value of the CV
    :param step_size: The size of the step for this CV
    :param num_bins: number of bins
    :param min_spring_constant: Minimal value of the spring constant to be used for this CV.
    :param max_spring_constant: Maximal value of the spring constant to be used for this CV.
    :param max_shift: Maximal value for the shift of this CV to compensate for free energy gradient.
    :param periodicity: periodicity of the CV
    :param units: the units of the CV
    :type cv_name: :class:`str`
    :type min_value: :class:`float`
    :type max_value: :class:`float`
    :type step_size: :class:`float`
    :type num_bins: :class:`int`
    :type periodicity: :class:`float`
    :type units: :class:`str`
    """
    self.name = cv_name
    self.min_value = min_value
    self.max_value = max_value
    self.min_spring_constant = min_spring_constant
    if max_spring_constant:
      self.max_spring_constant = max_spring_constant
    else:
      self.max_spring_constant = min_spring_constant
    if max_shift:
      self.max_shift = max_shift
    else:
      self.max_shift = 0.0
    self.step_size = step_size
    self.num_bins = num_bins
    self.bin_size = (self.max_value - self.min_value) / float(num_bins)
    self.periodicity = periodicity
    self.units = units
    # We extend the cv for the calculation of the PMF by step_size, so that windows
    # in the edge still have their free energy calculated properly
    if not self.periodicity:
      self.num_pads = 0
      self.wham_min_value = self.min_value - step_size
      self.wham_max_value = self.max_value + step_size
      self.wham_num_bins = num_bins + 2 * int(step_size / self.bin_size)
    else:
      self.num_pads = int(step_size / self.bin_size)
      self.wham_min_value = self.min_value
      self.wham_max_value = self.max_value
      self.wham_num_bins = num_bins
