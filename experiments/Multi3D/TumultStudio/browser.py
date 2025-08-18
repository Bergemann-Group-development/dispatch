from .widgets import EditableList
from .auxfuncs import *

class Browser(EditableList):
  def get_unique_label(self, label):
    labels = [o.label for o in self.objects]
    new_label = get_unique_name(label, labels, '_{:03d}')
    return new_label

  def get_objects_by_type(self, typ, exclude=None):
      objects = [o for o in self.objects if isinstance(o, typ)]
      if exclude:
          objects = [o for o in objects if not isinstance(o, exclude)]
      return objects
      
