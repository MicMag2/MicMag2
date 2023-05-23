import magnum.module as module
import magnum.magneto as magneto
from magnum.logger import logger
from magnum.mesh import VectorField, Field

class CurrentPath(module.Module):
  def __init__(self):
    super(CurrentPath, self).__init__()
    self.__contact1 = None
    self.__contact2 = None

  def calculates(self):
    return ["j", "j_potential", "R_contact"]

  def params(self):
    return ["sigma", "amr", "U_contact", "I_contact", "amr_dimension"]

  def on_param_update(self, id):
    # If material parameters change, we need re-setup self.__amr object.
    if id in ["sigma", "amr", "amr_dimension"]: self.__dirty = True
    # The following two parameters exclude each other, as contact_R is computed.
    if id == "U_contact": self.I_contact = None
    if id == "I_contact": self.U_contact = None

  def initialize(self, system):
    self.system = system
    self.sigma = Field(self.system.mesh); self.sigma.clear()
    self.amr = Field(self.system.mesh); self.amr.clear()

    # get contact structure from Bodies
    for i in range(len(system.world._World__bodies)):
      #print system.world._World__bodies[i].shape
      if system.world._World__bodies[i].id == "contact1":
        if (self.__contact1 != None): raise ValueError("CurrentPath module: 'contact1' defined more than once.")
        self.__contact1 = system.world._World__bodies[i].shape
      if system.world._World__bodies[i].id == "contact2":
        if (self.__contact2 != None): raise ValueError("CurrentPath module: 'contact2' defined more than once.")
        self.__contact2 = system.world._World__bodies[i].shape
    if (self.__contact1 == None) or (self.__contact2 == None) :
       raise ValueError("CurrentPath module: Not able to find two contact Bodies. Please create two instances of magnum.Body with id 'contact1' and 'contact2', respectively.")

    self.U_contact = 2.0
    self.I_contact = None
    self.amr_dimension = 2
    self.__dirty = True # self.__amr object still needs to be created and configured.

  def calculate(self, state, id):
    cache = state.cache

    def do_amr_setup(): # Re-setup AMR solver (self.__amr)
      # Defining cell types for AMR object: contact1_idxs, contact2_idxs, conductive_idxs index lists.
      contact1_idxs   = self.__contact1.getCellIndices(self.system.mesh)
      contact2_idxs   = self.__contact2.getCellIndices(self.system.mesh)
      conductive_idxs = [] # Put in indices of cells which have sigma!=0.0 here.
      for idx in range(self.system.mesh.total_nodes):
        if self.sigma.get(idx) != 0.0: conductive_idxs.append(idx)

      if contact1_idxs == []: raise ValueError("CurrentPath module: Cell range for contact 1 is empty. Try to slightly increase the size of the contact.")
      if contact2_idxs == []: raise ValueError("CurrentPath module: Cell range for contact 2 is empty. Try to slightly increase the size of the contact.")
      if contact1_idxs == contact2_idxs: raise ValueError("CurrentPath module: Contact pads are identical.")

      def list2intvector(lst): # Convert list of ints to std::vector<int> (called IntVector in Python) object. Needed for interfacing with C++ code.
        vec = magneto.IntVector()
        for i in lst: vec.push_back(i)
        return vec

      (nx, ny, nz), (dx, dy, dz) = self.system.mesh.num_nodes, self.system.mesh.delta
      self.__amr = magneto.AMR(nx, ny, nz, dx, dy, dz, list2intvector(conductive_idxs), list2intvector(contact1_idxs), list2intvector(contact2_idxs))
      self.__amr.setSigma(self.sigma)
      self.__amr.setAMR(self.amr, self.amr_dimension)

      # Mark as initialized and clear cache if neccessary.
      self.__dirty = False
      cache.amr_calculated = False
      if hasattr(cache, "amr_calculated"): delattr(cache, "amr_calculated")

    def do_amr(): # do_amr calculates the potential phi (stored interally in self.__amr) of the current path field j (for some internally chosen U_contact = const)
      if self.__dirty: do_amr_setup()
      if hasattr(cache, "amr_calculated"): return
      self.__amr.calculate(state.M)
      cache.amr_calculated = True

    def get_contact_voltage():
      # Get voltage 'U' applied to the contacts (in Volts).
      if state.U_contact: # Case 1: U is set explicitly
        U = state.U_contact(state.t) if callable(state.U_contact) else state.U_contact
      elif state.I_contact: # Case 2: I is set explicitly
        R = state.R_contact # This value is calculated from phi (see below) and cannot be set as a parameter.
        I = state.I_contact(state.t) if callable(state.I_contact) else state.I_contact
        U = I * R
      else: raise ValueError("CurrentPath: Need one of the parameters state.U_contact or state.I_contact.")
      return U

    if id == "j":
      if hasattr(cache, "amr_j"): return cache.amr_j
      j = cache.amr_j = VectorField(self.system.mesh)
      do_amr()
      self.__amr.get_j(j, get_contact_voltage())

      return j

    elif id == "j_potential":
      if hasattr(cache, "amr_j_potential"): return cache.amr_j_potential
      phi = cache.amr_j_potential = Field(self.system.mesh)
      do_amr()
      self.__amr.get_phi(phi, get_contact_voltage())
      return phi

    elif id == "R_contact":
      if hasattr(cache, "amr_R"): return cache.amr_R
      do_amr()
      R = cache.amr_R = self.__amr.get_resistance()
      return R

    else:
      raise KeyError("CurrentPath.calculate: Can't calculate %s", id)
